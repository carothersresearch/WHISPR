from whispr import *

"""
For user to change according to run
"""

# gene expression. 3ul rxns, 75% txtl
source_plate_type = '384PP_AQ_BP' 

sp_plasmids_file = 'source_plates/221004_arpae_plasmids.xlsx'
sp_plasmids = pd.read_excel(sp_plasmids_file, index_col = 0)
sp_plasmids = sp_plasmids[~sp_plasmids['Well'].isna()]

filename = '221004_arpae_genex.csv'
layout_genex_file = 'plate_layouts/'+filename 
layout_genex = pd.read_csv(layout_genex_file, index_col = 0, dtype = str)

mt_genex_file = 'mixing_tables/221004_arpae_plasmids.csv'
mt_genex = pd.read_csv(mt_genex_file, index_col = 0, dtype = str).fillna(0)

checkInputs(sp_plasmids, mt_genex,source_plate_type)

vol_table_df = generateVolumeTable(mt_genex, sp_plasmids, rxn_vol = 0.25*3, total_vol = 3, fill_with='water')

#specify rxn_vol (default = 2.5) and total_vol (default = 10) if you'd like to change the volume of each individual replicate or the total reaction volume


# biosynthesis. 25ul rxns, 2.5ul of diluted txtl

# buffer source plate
sp_buffers_file = 'source_plates/221004_arpae_buffers.xlsx'
sp_buffers= pd.read_excel(sp_buffers_file, index_col = 0)
sp_buffers = sp_buffers[~sp_buffers['Well'].isna()]

# HEPES reservoir 
sp_hepes = pd.DataFrame(columns=sp_buffers.columns)
sp_hepes = sp_hepes.append({'Label':'HEPES', 'Well':'A1,A2', 'Concentration':'', 'Volume':'2000,2000'}, ignore_index = True)

# gene expression source plate (diluted w/ Tris)
sp_genex = sp_from_layout(layout_genex, 60)

# group source plates
sp_types = ['384PP_AQ_BP','6RES_AQ_BP','384PP_AQ_BP'] # triple check this
sps = [sp_buffers,sp_hepes,sp_genex]

# get mixing table
mt_biosyn_file = 'mixing_tables/221004_arpae_buffers.csv'
mt_biosyn = pd.read_csv(mt_biosyn_file, index_col = 0, dtype = str).fillna(0)

# check formats
checkInputs(sps,mt_biosyn,sp_types)

# layouts for destination plate(s)
filename = '221004_arpae_biosyn_pos.csv'
layout_biosyn1_file = 'plate_layouts/'+filename 
layout_biosyn1 = pd.read_csv(layout_biosyn1_file, index_col = 0, dtype = str)
filename = '221004_arpae_biosyn_neg.csv'
layout_biosyn2_file = 'plate_layouts/'+filename 
layout_biosyn2 = pd.read_csv(layout_biosyn2_file, index_col = 0, dtype = str)
layouts = [layout_biosyn1,layout_biosyn2]

vol_table_df = generateVolumeTable(mt_biosyn, sps, rxn_vol = 22.5, total_vol = 25, fill_with='HEPES')

protocol_biosyn_df = writeProtocol(sp_types, vol_table_df, layouts,sps, update_source_vol=False)

# protocol_biosyn_df.to_csv('protocols/221004_arpae_biosyn.csv',index = False)