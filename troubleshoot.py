from whispr import *

"""
For user to change according to run
"""

#CSV file with source plate info. Columns are label, well, concentration, volume
source_plate = 'Experiments/CRISPRa_degradation/230221_sp.xlsx'
source_plate_df = pd.read_excel(source_plate, index_col = 0, convert_float = False)
source_plate_df = source_plate_df[~source_plate_df['Well'].isna()]

filename = '230217_pl.csv'
# CSV file with output plate layout: columns labeled 1-12 and rows labeled A-H for 96 well plate
output_layout = 'Experiments/CRISPRa_degradation/'+filename 
output_plate = pd.read_csv(output_layout, index_col = 0, dtype = str)

#CSV file with each reaction as rows and the volume of each input (columns) added to reach final reaction volume
filename = '230217_mt.csv'
mixing_table = 'Experiments/CRISPRa_degradation/'+filename 
mixing_table_df = pd.read_csv(mixing_table, index_col = 0, dtype = str).fillna(0)

source_plate_type = '384PP_AQ_BP' 

checkInputs(source_plate_df, mixing_table_df,source_plate_type)

vol_table_df = generateVolumeTable(mixing_table_df, source_plate_df, rxn_vol = 50, total_vol = 50, fill_with=False)

#specify rxn_vol (default = 2.5) and total_vol (default = 10) if you'd like to change the volume of each individual replicate or the total reaction volume
filename = '230221_ep.csv'
output_df = writeProtocol(source_plate_type, vol_table_df, output_plate,source_plate_df)[0]
output_df.to_csv('Experiments/CRISPRa_degradation/'+filename,index = False)