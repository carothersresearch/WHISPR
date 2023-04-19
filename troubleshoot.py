from whispr import *
import os 

"""
For user to change according to run
"""
# gene expression. 2.5ul rxns, 75% txtl
rxn_vol = 3
source_plate_type = '384PP_AQ_BP' 

os.getcwd()
folder = os.getcwd() + '/Experiments/230419_mod1-3/'
sp_plasmids_file = folder + 'plasmids_sp_230419.xlsx'
sp_plasmids = pd.read_excel(sp_plasmids_file, index_col = 0, engine='openpyxl')
sp_plasmids = sp_plasmids[~sp_plasmids['Well'].isna()]

layout_genex_file = folder + 'genex-pl-384.csv'
layout_genex = pd.read_csv(layout_genex_file, index_col = 0, dtype = str)

mt_genex_file = folder + 'genex_mt_230419.xlsx'
mt_genex = pd.read_excel(mt_genex_file, index_col = 0, dtype = str, engine = 'openpyxl').fillna(0)

# biosynthesis. 25ul rxns, 2.5ul of diluted txtl

# buffer source plate
sp_buffers_file = folder + 'buffers_sp.xlsx'
sp_buffers= pd.read_excel(sp_buffers_file, index_col = 0, engine='openpyxl')
sp_buffers = sp_buffers[~sp_buffers['Well'].isna()]

# HEPES reservoir (actually same plate as Tris)
sp_hepes = pd.DataFrame(columns=sp_buffers.columns)
sp_hepes = sp_hepes.append({'Label':'HEPES', 'Well':'A2,A3', 'Concentration':'', 'Volume':'2000,2000'}, ignore_index = True)

# gene expression source plate (diluted w/ Tris)
sp_genex = sp_from_layout(layout_genex, 60)

# group source plates
sp_types = ['384PP_AQ_BP','6RES_AQ_BP2','384PP_AQ_BP'] # triple check this
sps = [sp_buffers,sp_hepes,sp_genex]

# get mixing table
mt_biosyn_file = folder + 'buffers-mt-noHT.xlsx'
mt_biosyn = pd.read_excel(mt_biosyn_file, index_col = 0, dtype = str, engine = 'openpyxl').fillna(0)

# check formats
checkInputs(sps,mt_biosyn,sp_types)

# layouts for destination plate(s)
layout_biosyn1_file = folder + 'biosyn_pl.csv' 
layout_biosyn1 = pd.read_csv(layout_biosyn1_file, index_col = 0, dtype = str)
# filename = '221018_arpae_biosyn_neg.csv'
# layout_biosyn2_file = 'plate_layouts/'+filename 
# layout_biosyn2 = pd.read_csv(layout_biosyn2_file, index_col = 0, dtype = str)
# layouts = [layout_biosyn1,layout_biosyn2]

vol_table_df = generateVolumeTable(mt_biosyn, sps, rxn_vol = 22.5, total_vol = 25, fill_with='HEPES')

protocol_biosyn_dfs = writeProtocol(sp_types, vol_table_df, layout_biosyn1,sps, update_source_vol= folder +'combined_sp_updated.xlsx')

protocol_biosyn_dfs[0].to_csv(folder + 'biosyn_protocol_384.csv',index = False)
protocol_biosyn_dfs[1].to_csv(folder + 'biosyn_protocol_6RES.csv',index = False)