#plate_type, vol_table, source_plate_layout, output_layout

#CSV file with source plate info. Columns are label, well, concentration, volume
source_plate = 'C:/users/rycar/Desktop/UW/Carothers/Python/ECHO_source_plate.xlsx'
source_plate_df = pd.read_excel(source_plate, index_col = 0)

# CSV file with output plate layout: columns labeled 1-12 and rows labeled A-H for 96 well plate
output_layout = 'C:/users/rycar/Desktop/UW/Carothers/Plate_Reader/plate_layout_test.csv' 

#CSV file with each reaction as rows and the volume of each input (columns) added to reach final reaction volume
mixing_table = 'C:/users/rycar/Desktop/UW/Carothers/Python/mixing_table_conc_test.csv' 
mixing_table_df = pd.read_csv(mixing_table, index_col = 0, dtype = str)

source_plate_type = '384PP_AQ_BP' 

vol_table = generateVolumeTable(mixing_table_df, source_plate_df)

'''
-------------------------------------------------------------
'''

# check source plate type and set volume range
if 'LDV' in plate_type:
    vol_range = 9.5
elif '384PP' in plate_type:
    vol_range = 45

# keep track of volume used for each component
vol_used = {}
for k in list(vol_table.columns[1:]):
    vol_used[k] = 0


# reads plate layout and assigns wells to each reaction (rxn_loc; dict)
plate = pd.read_csv(output_layout, index_col = 0, dtype = str) 
labels = pd.unique(np.concatenate(plate.values))
rxn_loc = {}
for l in labels:
    if type(l) is str:
        index = plate[plate.isin([l])].stack().index
        rxn_loc[l] = []
        for i in index:
            rxn_loc[l].append(str(i[0]) + i[1])



# create output dataframe
output_df = pd.DataFrame(columns = {'Source Plate Name', 'Source Plate Type', 'Source Well', 
                                 'Destination Plate Name', 'Destination Well', 'Transfer Volume'})


'''
for each reaction in the plate layout:
    if this reaction label matches one in the mixing table:
        for each well of the reaction in the plate layout:
            remove "total" volume column
            for each component in the reaction:
                find volume to transfer of component
                if any volume is added:
                    if the volume used of this component is less than the volume range:
                        write row in protocol using first source well for component
                    if the volume used of this component is greater than the volume range:
                        write row in protocol using second source well for component
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
                ## separate if there is > 1 well in source plate
                    source_well = source_plate_df.loc[component]['Well']
                    if type(source_well) == str:
                    ## only one source well
                        if vol_used[component] + transfer_vol < vol_range:
                            row = {'Source Plate Name':'Source[1]', 'Source Plate Type': plate_type, 'Source Well': source_well,
                                'Destination Plate Name':'Destination[1]', 'Destination Well': well, 'Transfer Volume': transfer_vol*1000}
                            vol_used[component] = vol_used[component] + transfer_vol

                        else:
                            raise NameError('Need more volume of ' +component+ ' to complete reaction. Add another well to source plate.')

                    else:
                    ## more than one source well, will be a panda series
                    ## convert to list
                        source_well = list(source_well)
                        ## use first well unless the well is empty, then use second well

                        if vol_used[component] + transfer_vol < vol_range:
                            row = {'Source Plate Name':'Source[1]', 'Source Plate Type': plate_type, 'Source Well': source_well[0],
                                'Destination Plate Name':'Destination[1]', 'Destination Well': well, 'Transfer Volume': transfer_vol*1000}

                        elif vol_used[component] + transfer_vol >= vol_range:
                            # do i still need this?
                            source_well = source_well[1:]
                            if len(source_well) == 0:
                                raise NameError('Need more volume of ' +component+ ' to complete reaction. Add another well to source plate.')
                            row = {'Source Plate Name':'Source[1]', 'Source Plate Type': plate_type, 'Source Well': source_well[0],
                                'Destination Plate Name':'Destination[1]', 'Destination Well': well, 'Transfer Volume': transfer_vol*1000}

                            vol_used[component] = 0

                        vol_used[component] = vol_used[component] + transfer_vol


                    output_df = output_df.append(row, ignore_index = True)
                    

output_df = output_df[['Source Plate Name', 'Source Plate Type', 'Source Well', 
                        'Destination Plate Name', 'Destination Well', 'Transfer Volume']]





print('\n')
for component in vol_used:
    if 'LDV' in plate_type:
        print('Load at least ', np.round(4.5+vol_used[component],2), 'ul and maximum 14 ul of ', component)
    elif '384PP' in plate_type:
        print('Load at least ', np.round(20+vol_used[component],2), 'ul and maximum 65 ul of ', component)


