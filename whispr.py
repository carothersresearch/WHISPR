import pandas as pd
import numpy as np

#This function rounds to the nearest .025 as the ECHO only adds in .025ul increments.
def myround(x, prec=3, base=.025):
    return round(base * round(float(x)/base),prec)

def checkInputs(source_plate, mixing_table_df, plate_type = '384PP_AQ_BP'):

    '''

    Checks that all of the inputs are in the correct volume range for source plate type
    
    Parameters:
    -----------
    source plate: dataframe with source plate info. Columns are label, well, concentration, volume
    mixing table: dataframe with each reaction as rows and the volume of each input (columns) added to reach final reaction volume
    plate type (opt): (default = 384PP_AQ_BP) specifiy source plate type
    
    Output:
    --------
    None, error if outside of working volume range


    '''

    # check source plate type and set volume range

    if not source_plate['Well'].is_unique:
        raise NameError('Wells in the source plate are not unique!')

    present = [m in source_plate.index for m in mixing_table_df.columns]
    if not np.all(present):
        raise NameError('Source plate does not contain some plasmids in the mixing table: '+ str(mixing_table_df.columns[not present]))

    if 'LDV' in plate_type:
        vol_min = 4.5
        vol_max = 14
        vol_range = vol_max - vol_min

    elif '384PP' in plate_type:
        vol_min = 19
        vol_max = 65
        vol_range = vol_max - vol_min

        
    vol = []
    for component in list(source_plate['Volume']):
        if type(component) == str:
            c = component.split(',')
            for i in c:
                vol.append(float(i))
        else:
            vol.append(component)
            
    if any([v > vol_max for v in vol]):
        raise NameError('Volumes of source plate are above working volume range.')
    if any([v < vol_min for v in vol]):
        raise NameError('Volumes of source plate are below working volume range.')






def generateVolumeTable(mixing_table_df, source_plate_df, rxn_vol = 2.5, total_vol = 10):
    '''

    Converts concentrations to volumes for reaction mixing table, raises error if volume is above max (rxn_vol, default = 2.5ul)
    
    Parameters:
    -----------
    source plate: dataframe with source plate info. Columns are label, well, concentration, volume
    mixing table: dataframe with each reaction as rows and the volume of each input (columns) added to reach final reaction volume
    
    Output:
    --------
    Dataframe with volumes to add of each component for each reaction


    '''



    vol_table = []
    vol_table_df = pd.DataFrame(columns = ['Label'] + list(source_plate_df['Label']))

    for row in mixing_table_df.index:

        vol_table_df = vol_table_df.append({'Label': str(row)}, ignore_index = True)    
        vol = 0

        for column in mixing_table_df.columns:
            conc_to_add = float(mixing_table_df.loc[row][column])
            label_indx = 0
            conc_of_source = source_plate_df.loc[column]['Concentration']
            if type(conc_of_source) != np.float64:
                conc_of_source = conc_of_source.sort_values(ascending = False)[label_indx]
                
            #conc_of_source = source_plate_df[source_plate_df['Label'] == column]['Concentration'].sort_values(ascending = False)[label_indx]
            vol_to_add = myround(total_vol*conc_to_add/conc_of_source)
            while conc_to_add > 0 and vol_to_add == 0:
                if type(source_plate_df.loc[column]['Concentration']) == np.float64:
                    raise NameError('Mate you need a more dilute stock of '+column)
                else:
                    label_indx += 1
                    if label_indx >= len(source_plate_df.loc[column]['Concentration']):
                        raise NameError('Mate you need a more dilute stock of '+column)
                    else:
                        conc_of_source = source_plate_df.loc[column]['Concentration'].sort_values(ascending = False)[label_indx]
                vol_to_add = myround(total_vol*conc_to_add/conc_of_source)

            label = source_plate_df[source_plate_df['Concentration'] == conc_of_source].loc[column]['Label']
           # label = source_plate_df[source_plate_df['Concentration'] == conc_of_source]
            vol_table_df.loc[vol_table_df['Label'] == row,label] = vol_to_add

            vol+=vol_to_add

        if round(vol,3) > rxn_vol:
            print(vol_table_df)
            raise NameError('Volume of '+ row+ ' exceeds '+str(rxn_vol) + 'ul. Total volume is '+ str(round(vol,3))+' Please change volumes and try again.')
        else:

            vol_table_df.loc[vol_table_df['Label'] == row,'Water'] = myround(rxn_vol - vol)
            vol_table_df.loc[vol_table_df['Label'] == label,column] = vol_to_add



    return vol_table_df



def writeProtocol(plate_type, vol_table, output_layout,source_plate_df, update_source_vol = None):
        '''
        Writes protocol for use with ECHO plate reader

        Parameters:
        ----------
        - plate_type: Source plate calibration (str)
        - mixing_table: path to csv file with reaction volumes (str)
        - input_layout: links inputs to source well (dict)
        - output_layout: path to csv with desired plate layout for 96 well plate (str)
        
        
        Returns:
        --------
        - dataframe of ECHO protocol


        Reference slides here for more information: 
        https://docs.google.com/presentation/d/1VzEFFiyCCfI-mrfSQGjb41TOOcz61sTlfk1WMnb-7BI/edit#slide=id.gf541592c34_0_28
        '''

        

        

        # check source plate type and set volume range
        if 'LDV' in plate_type:
            vol_min = 4.5
            vol_max = 14
            vol_range = vol_max - vol_min

        elif '384PP' in plate_type:
            vol_min = 20
            vol_max = 65
            vol_range = vol_max - vol_min

        # keep track of volume used for each component
        vol_used = {}
        well_list = ''
        for k in list(source_plate_df['Well']):
            if type(k) is str:
                well_list += k.replace(' ','') + ','
        well_list = well_list.split(',')
        for k in well_list:
            vol_used[k] = 0


        # reads plate layout and assigns wells to each reaction (rxn_loc; dict)
        labels = pd.unique(np.concatenate(output_layout.values))
        rxn_loc = {}
        for l in labels:
            if type(l) is str:
                index = output_layout[output_layout.isin([l])].stack().index
                rxn_loc[l] = []
                for i in index:
                    rxn_loc[l].append(str(i[0]) + str(i[1]))



        # create output dataframe
        output_df = pd.DataFrame(columns = {'Source Plate Name', 'Source Plate Type', 'Source Well', 
                                        'Destination Plate Name', 'Destination Well', 'Transfer Volume'})


        '''
        for each reaction in the plate layout:
            if this reaction label matches one in the mixing table:
                for each well of the reaction in the plate layout:
                    for each component in the reaction:
                        find volume to transfer of component
                        if any volume is added:
                            if there is one source well:
                                if the volume used of this component is less than the volume range:
                                    append row
                            if there is more than one source well:
                                if the volume used of this component is less than the volume range:
                                    append row
                                else:
                                    remove first source well
                                    append row
                            update volume used
                            append output dataframe                    

        ''' 


        # subtract from volume in source plate file 
        well_vols = {}
        for component in list(source_plate_df['Label']):
            well_vols[component] = str(list(source_plate_df[source_plate_df['Label'] == component]['Volume'])[0]).split(',')

        rxn_keys = list(rxn_loc.keys())
        for rxn in rxn_keys:
            if rxn in list(vol_table['Label']): 
                for well in rxn_loc[rxn]: 
                    vols = vol_table[vol_table['Label'] == rxn]
                    for component in vols.columns[1:]: 
                        transfer_vol = float(vols[component])      
                        if transfer_vol > 0:
                        ## separate if there is > 1 well in source plate
                            source_well = list(source_plate_df[source_plate_df['Label'] == component]['Well'])
                            if not type(source_well[0]) == list: source_well = source_well[0].split(',')

                            ## use first well unless the well is empty, then use second well

                            # check if volume used leaves volume below minimum
                            if vol_used[source_well[0]] + transfer_vol >= float(well_vols[component][0]) - vol_min:
                                print(source_well)
                                if len(source_well) <= 1:
                                    raise NameError('Need more volume of ' +component+ ' to complete reaction ' + rxn + '. Add another well to source plate.')

                                source_plate_df['Well'][component] = ','.join(source_well[1:]).replace(' ','')
                                well_vols[component] = well_vols[component][1:]
                                
                            row = {'Source Plate Name':'Source[1]', 'Source Plate Type': plate_type, 'Source Well': source_well[0],
                                'Destination Plate Name':'Destination[1]', 'Destination Well': well, 'Transfer Volume': transfer_vol*1000}

                            vol_used[source_well[0]] = vol_used[source_well[0]] + transfer_vol

                            output_df = output_df.append(row, ignore_index = True)


        output_df = output_df[['Source Plate Name', 'Source Plate Type', 'Source Well', 
                                'Destination Plate Name', 'Destination Well', 'Transfer Volume']]


        for v in vol_used:
            vol_used[v] = myround(vol_used[v])
        vol_used = {well:vol for well,vol in vol_used.items() if vol!=0}
        print('Volumes used from each well for this protocol:')
        print(vol_used)

        if update_source_vol:
            ignore_water = True # at least for now. might get to it later
            if ignore_water:
                for k,v in vol_used.items():
                    if 'A' not in k:
                        row = 1 + list(source_plate_df['Well'][1:].values).index(k)
                        current_vol = source_plate_df.iloc[row,3]
                        new_vol = current_vol - v
                        source_plate_df.iloc[row,4] = new_vol

            col_names = list(source_plate_df.columns)
            col_names[4] = 'New Volume'
            source_plate_df.columns = col_names
            source_plate_df.to_excel(update_source_vol)

        return output_df