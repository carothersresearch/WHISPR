from tokenize import Name
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
    if type(source_plate) is not list: source_plate = [source_plate]
    if type(plate_type) is not list: plate_type = [plate_type]
    if len(source_plate) is not len(plate_type): 
        raise NameError('Mismatch in number of plates and plate types')
    else:
        for sp,spt,k in zip(source_plate,plate_type,range(len(source_plate))):

        # check source plate type 

            if not sp['Well'].is_unique:
                raise NameError('Wells in the source plate are not unique!')

            if 'LDV' in spt:
                vol_min = 4.5
                vol_max = 14

            elif '384PP' in spt:
                vol_min = 18
                vol_max = 65

            elif '6RES' in spt:
                vol_min = 250
                vol_max = 2800

            vol = []
            for component in list(sp['Volume']):
                if type(component) == str:
                    c = component.split(',')
                    for i in c:
                        vol.append(float(i))
                else:
                    vol.append(component)
                    
            if any([v > vol_max for v in vol]):
                raise NameError('Volumes of source plate '+str(k)+' are above working volume range.')
            if any([v < vol_min for v in vol]):
                raise NameError('Volumes of source plate '+str(k)+' are below working volume range.')
        
        if not np.all([m in combine_sps(source_plate)['Label'] for m in mixing_table_df.columns]):
            raise NameError('Source plate does not contain some items in the mixing table')

        for sp in source_plate:
            if sp.index.name != 'Item':
                sp.index = sp['Label']
                sp.index.name = 'Item'

def generateVolumeTable(mixing_table_df, source_plate_df, rxn_vol = 2.5, total_vol = 10, min_well_vol = 18, fill_with = 'Water', multiRpW = False):
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
    if type(source_plate_df) is list:
        source_plate_df = combine_sps(source_plate_df)

    running_source_plate = source_plate_df.copy(deep = True)

    vol_table = []
    vol_table_df = pd.DataFrame(columns = ['Label'] + list(running_source_plate['Label']))

    for row in mixing_table_df.index:

        vol_table_df = vol_table_df.append({'Label': str(row)}, ignore_index = True)    
        vol = 0

        for column in mixing_table_df.columns:
            conc_to_add = float(mixing_table_df.loc[row][column])
            label_indx = 0
            conc_of_source = running_source_plate.loc[column]['Concentration']
            if (type(conc_of_source) != np.float64) and (type(conc_of_source) != float):
                while running_source_plate.loc[column].sort_values(ascending = False, by = 'Concentration')['Volume'][label_indx] <= min_well_vol:
                    label_indx += 1
                else: 
                    conc_of_source = conc_of_source.sort_values(ascending = False)[label_indx]
                    
            #conc_of_source = source_plate_df[source_plate_df['Label'] == column]['Concentration'].sort_values(ascending = False)[label_indx]
            vol_to_add = myround(total_vol*conc_to_add/conc_of_source)
            # this may round to zero, so need to check for more dilute wells 
            while conc_to_add > 0 and vol_to_add == 0:
                if type(running_source_plate.loc[column]['Concentration']) == np.float64:
                    raise NameError('Mate you need a more dilute stock of '+column)
                else:
                    label_indx += 1
                    if label_indx >= len(running_source_plate.loc[column]['Concentration']):
                        raise NameError('Mate you need a more dilute stock of '+column)
                    else:
                        conc_of_source = running_source_plate.loc[column]['Concentration'].sort_values(ascending = False)[label_indx]
                vol_to_add = myround(total_vol*conc_to_add/conc_of_source)

            label = running_source_plate[running_source_plate['Concentration'] == conc_of_source].loc[column]['Label']
           # label = source_plate_df[source_plate_df['Concentration'] == conc_of_source]
            vol_table_df.loc[vol_table_df['Label'] == row,label] = vol_to_add
            
            running_source_plate=running_source_plate.reset_index(level=0).set_index('Label')
            if type(running_source_plate.loc[label, 'Volume']) == str: 
                current_vols = np.array(list(map(float,running_source_plate.loc[label, 'Volume'].split(','))))
                current_vols[np.where(current_vols>min_well_vol)[0][0]] = current_vols[np.where(current_vols>min_well_vol)[0][0]] - vol_to_add
                running_source_plate.loc[label, 'Volume'] = ','.join(list(map(str,current_vols)))
            else:
                running_source_plate.loc[label, 'Volume'] = running_source_plate.loc[label, 'Volume'] - vol_to_add
            running_source_plate=running_source_plate.reset_index(level=0).set_index('Item') 
            vol+=vol_to_add

        if round(vol,3) > rxn_vol:
            print(vol_table_df.loc[vol_table_df['Label'] == row].to_markdown())
            raise NameError('Volume of '+ row+ ' exceeds '+str(rxn_vol) + 'ul. Total volume is '+ str(round(vol,3))+' Please change volumes and try again.')
        else:
            if fill_with:
                vol_table_df.loc[vol_table_df['Label'] == row,fill_with] = myround(rxn_vol - vol)
                
                running_source_plate=running_source_plate.reset_index(level=0).set_index('Label')
                if type(running_source_plate.loc[fill_with, 'Volume']) == str: 
                    current_vols = np.array(list(map(float,running_source_plate.loc[fill_with, 'Volume'].split(','))))
                    current_vols[np.where(current_vols>min_well_vol)[0][0]] = current_vols[np.where(current_vols>min_well_vol)[0][0]] - myround(rxn_vol - vol)
                    running_source_plate.loc[fill_with, 'Volume'] = ','.join(list(map(str,current_vols)))
                else:
                    running_source_plate.loc[fill_with, 'Volume'] = running_source_plate.loc[fill_with, 'Volume'] - myround(rxn_vol - vol)
                running_source_plate=running_source_plate.reset_index(level=0).set_index('Item')
    
    if multiRpW:
        for i in range(len(vol_table_df)):
            try: 
                nrxn = int(vol_table_df['Label'][i].split('_')[0])
            except:
                nrxn = 1
            vol_table_df.iloc[i,1:] = vol_table_df.iloc[i,1:]*nrxn

    return vol_table_df

def sp_from_layout(layout, volumes):
    sp = pd.DataFrame()

    labels = pd.unique(np.concatenate(layout.values))
    exps = {}
    xlbl = []
    for l in labels:
        if type(l) is str:
            xlbl.append(l)
            index = layout[layout.isin([l])].stack().index
            exps[l] = []
            for i in index:
                exps[l].append(str(i[0]) + i[1])

    if (type(volumes) is int) or (type(volumes) is float):
        for label, wells in exps.items():
            sp = sp.append({'Label':label, 'Well':','.join(wells), 'Concentration':100, 'Volume':','.join([str(volumes)]*len(wells))}, ignore_index = True)
    else:
        raise NameError('volumes type not supported yet')

    sp.index = sp['Label']
    sp.index.name = 'Item'
    return sp

def combine_sps(list_source_plates):
    sp = pd.concat(list_source_plates, keys=range(1,len(list_source_plates)+1)).reset_index(level=0).rename(columns={'level_0': 'PlateID'})
    if not sp['Label'].is_unique:
        raise NameError('Items in the source plates are not unique!')
    sp.index.name = 'Item'
    return sp

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

        if type(source_plate_df) is list:
            if len(source_plate_df) is not len(plate_type): 
                raise NameError('Mismatch in number of plates and plate types')
        else:
            source_plate_df = [source_plate_df]
            plate_type = [plate_type]

        source_plate_df = combine_sps(source_plate_df)


        running_source_plate = source_plate_df.copy(deep = True)


        # check source plate type and set volume range
        def get_vol_min(plate_type):
            if 'LDV' in plate_type:
                vol_min = 4.5

            elif '384PP' in plate_type:
                vol_min = 18 # 20

            elif '6RES' in plate_type:
                vol_min = 250

            return vol_min

        # keep track of volume used for each component
        vol_used = {}
        well_list = ''
        for p,k in zip(running_source_plate['PlateID'],list(running_source_plate['Well'])):
            p = str(p)
            if type(k) is str:
                well_list += p+'_'+k.replace(' ','').replace(',',','+p+'_') + ','
        well_list = well_list[:-1].split(',')

        vol_list = ''
        for k in list(running_source_plate['Volume']):
            if (type(k) is int) or (type(k) is float):
                vol_list += str(k) + ','
            if type(k) is str:
                vol_list += k.replace(' ','') + ','
        vol_list = vol_list[:-1].split(',')

        for k in well_list:
            vol_used[k] = 0

        vol_left = dict(zip(well_list, map(float,vol_list)))

        # subtract from volume in source plate file 
        well_vols = {}
        for component in list(running_source_plate['Label']):
            well_vols[component] = str(list(running_source_plate[running_source_plate['Label'] == component]['Volume'])[0]).split(',')

        # create output dataframe
        output_df = pd.DataFrame(columns = {'Source Plate Name', 'Source Plate Type', 'Source Well', 
                                        'Destination Plate Name', 'Destination Well', 'Transfer Volume'})
                                        

        if type(output_layout) is not list:
            output_layout = [output_layout]

        for layout_id, layout in enumerate(output_layout):    
            # reads plate layout and assigns wells to each reaction (rxn_loc; dict)
            labels = pd.unique(np.concatenate(layout.values))
            rxn_loc = {}
            for l in labels:
                if type(l) is str:
                    index = layout[layout.isin([l])].stack().index
                    rxn_loc[l] = []
                    for i in index:
                        rxn_loc[l].append(str(i[0]) + str(i[1]))



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

            rxn_keys = list(rxn_loc.keys())
            for rxn in rxn_keys:
                if rxn in list(vol_table['Label']): 
                    for well in rxn_loc[rxn]: 
                        vols = vol_table[vol_table['Label'] == rxn]
                        for component in vols.columns[1:]: 
                            transfer_vol = float(vols[component])      
                            if transfer_vol > 0:
                                row = running_source_plate[running_source_plate['Label'] == component]
                            ## separate if there is > 1 well in source plate
                                source_well = list(row['Well'])
                                if not type(source_well[0]) == list: source_well = source_well[0].split(',')

                                ## use first well unless the well is empty, then use second well
                                index = 0
                                
                                # check if volume used leaves volume below minimum
                                if vol_used[str(row['PlateID'][0])+'_'+source_well[0]] + transfer_vol >= float(well_vols[component][0]) - get_vol_min(plate_type[row['PlateID'][0]-1]):
                                    # print(vol_used[source_well[0]], transfer_vol)
                                    # print(source_well)
                                    if len(source_well) <= 1:
                                        raise NameError('Need more volume of ' +component+ ' to complete reaction ' + rxn + '. Add another well to source plate.')

                                    running_source_plate.loc[running_source_plate['Label'] == component,'Well'] = ','.join(source_well[1:]).replace(' ','')
                                    well_vols[component] = well_vols[component][1:]
                                    index = 1
    
                                echo_command = {'Source Plate Name':'Source['+str(row['PlateID'][0])+']', 'Source Plate Type': plate_type[row['PlateID'][0]-1], 'Source Well': source_well[index],
                                    'Destination Plate Name':'Destination['+str(layout_id+1)+']', 'Destination Well': well, 'Transfer Volume': transfer_vol*1000}

                                vol_used[str(row['PlateID'][0])+'_'+source_well[index]] = vol_used[str(row['PlateID'][0])+'_'+source_well[index]] + transfer_vol

                                output_df = output_df.append(echo_command, ignore_index = True)


        output_df = output_df[['Source Plate Name', 'Source Plate Type', 'Source Well', 
                                'Destination Plate Name', 'Destination Well', 'Transfer Volume']]


        for v in vol_used:
            vol_used[v] = myround(vol_used[v])
            vol_left[v] -= vol_used[v]

        vol_used = {well:vol for well,vol in vol_used.items() if vol!=0}
        print('Volumes used from each well for this protocol:')
        print(vol_used)

        new_vol_list = ['']*len(source_plate_df)
        for k,v in vol_left.items():
            row = (source_plate_df['PlateID']==int(k[0])*source_plate_df['Well'].str.contains(k[2:])).argmax()
            new_vol_list[row] = new_vol_list[row] + str(v) + ','

        source_plate_df['New Volume'] = new_vol_list
        if update_source_vol:
            source_plate_df.to_excel(update_source_vol)

        output_df.sort_values(['Destination Plate Name', 'Source Plate Name'])
        outputs = []
        for pt in np.unique(plate_type):
            outputs.append(output_df[output_df['Source Plate Type'] == pt])
        
        return outputs