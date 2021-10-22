from whispr import *

#CSV file with source plate info. Columns are label, well, concentration, volume
source_plate = 'source_plates/ECHO_source_plate.xlsx'
source_plate_df = pd.read_excel(source_plate, index_col = 0)

filename = '211021-titrations-test1.csv'
# CSV file with output plate layout: columns labeled 1-12 and rows labeled A-H for 96 well plate
output_layout = 'plate_layouts/'+filename 
output_plate = pd.read_csv(output_layout, index_col = 0, dtype = str)

#CSV file with each reaction as rows and the volume of each input (columns) added to reach final reaction volume
mixing_table = 'mixing_tables/'+filename 
mixing_table_df = pd.read_csv(mixing_table, index_col = 0, dtype = str).fillna(0)

source_plate_type = '384PP_AQ_BP' 

checkInputs(source_plate_df, mixing_table_df,source_plate_type)

vol_table_df = generateVolumeTable(mixing_table_df, source_plate_df)

output_df = writeProtocol(source_plate_type, vol_table_df, source_plate_df, output_plate,source_plate_df)
output_df.to_csv('protocols/'+filename,index = False)