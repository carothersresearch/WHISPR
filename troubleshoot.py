from whispr import *

"""
For user to change according to run
"""

# gene expression. 3ul rxns, 75% txtl
source_plate_type = '384PP_AQ_BP' 

sp_plasmids_file = 'source_plates/221018_arpae_plasmids.xlsx'
sp_plasmids = pd.read_excel(sp_plasmids_file, index_col = 0)
sp_plasmids = sp_plasmids[~sp_plasmids['Well'].isna()]

filename = '221018_arpae_genex.csv'
layout_genex_file = 'plate_layouts/'+filename 
layout_genex = pd.read_csv(layout_genex_file, index_col = 0, dtype = str)

mt_genex_file = 'mixing_tables/221018_arpae_plasmids.csv'
mt_genex = pd.read_csv(mt_genex_file, index_col = 0, dtype = str).fillna(0)

checkInputs(sp_plasmids, mt_genex,source_plate_type)

vol_table_df = generateVolumeTable(mt_genex, sp_plasmids, rxn_vol = 0.25*3, total_vol = 3, fill_with='water', multiRpW = True)

#specify rxn_vol (default = 2.5) and total_vol (default = 10) if you'd like to change the volume of each individual replicate or the total reaction volume
protocol_genex_df = writeProtocol(source_plate_type, vol_table_df, layout_genex,sp_plasmids, update_source_vol='source_plates/221018_arpae_plasmids.xlsx')
