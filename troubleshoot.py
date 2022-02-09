from whispr import *

"""
For user to change according to run
"""

#CSV file with source plate info. Columns are label, well, concentration, volume
source_plate = 'source_plates/220208_source_plate.xlsx'
source_plate_df = pd.read_excel(source_plate, index_col = 0)
source_plate_df = source_plate_df[~source_plate_df['Well'].isna()]

filename = '220209_split_scRNAs_PL.csv'
# CSV file with output plate layout: columns labeled 1-12 and rows labeled A-H for 96 well plate
output_layout = 'plate_layouts/'+filename 
output_plate = pd.read_csv(output_layout, index_col = 0, dtype = str)

#CSV file with each reaction as rows and the volume of each input (columns) added to reach final reaction volume
filename = '220209_split_scRNAs_MT.csv'
mixing_table = 'mixing_tables/'+filename 
mixing_table_df = pd.read_csv(mixing_table, index_col = 0, dtype = str).fillna(0)

source_plate_type = '384PP_AQ_BP' 

checkInputs(source_plate_df, mixing_table_df,source_plate_type)

vol_table_df = generateVolumeTable(mixing_table_df, source_plate_df)

#specify rxn_vol (default = 2.5) and total_vol (default = 10) if you'd like to change the volume of each individual replicate or the total reaction volume
filename = '220202_split_scRNAs_EP.csv'
output_df = writeProtocol(source_plate_type, vol_table_df, output_plate,source_plate_df, update_source_vol = source_plate)
output_df.to_csv('protocols/'+filename,index = False)