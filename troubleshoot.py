from whispr import *

"""
For user to change according to run
"""


file_name = "C:/Users/rycar/Documents/GitHub/WHISPR/inputs/211105_BL-CRISPRa.xlsx"
#CSV file with source plate info. Columns are label, well, concentration, volume
source_plate_df = pd.read_excel(file_name, sheet_name = 'ECHO_source_plate', index_col = 0)

# CSV file with output plate layout: columns labeled 1-12 and rows labeled A-H for 96 well plate
output_plate = pd.read_excel(file_name, sheet_name = "plate_layout", index_col = 0, dtype = str)

#CSV file with each reaction as rows and the volume of each input (columns) added to reach final reaction volume
mixing_table_df = pd.read_excel(file_name, sheet_name = "mixing_table", index_col = 0, dtype = str).fillna(0)

source_plate_type = '384PP_AQ_BP' 

checkInputs(source_plate_df, mixing_table_df,source_plate_type)

#specify rxn_vol (default = 2.5) and total_vol (default = 10) if you'd like to change the volume of each individual replicate or the total reaction volume
vol_table_df = generateVolumeTable(mixing_table_df, source_plate_df,rxn_vol = 2.5)

output_df = writeProtocol(source_plate_type, vol_table_df, source_plate_df, output_plate,source_plate_df)
filename = '211105-BL-CRISPRa_1.csv'
output_df.to_csv('protocols/'+filename,index = False)