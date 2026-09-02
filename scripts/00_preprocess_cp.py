"""Preprocess iDC1003 for joint-model integration:
1. Add transport reactions for missing critical metabolites (formate, GlcNAc, fucose L, urea)
2. Build CP-specific id_mapping (CA mapping reused + new additions + D-lactate/NH3 specials)
3. Save modified model
"""
import cobra
import pandas as pd
from cobra import Metabolite, Reaction

print('Loading iDC1003...')
cp = cobra.io.read_sbml_model('models/cp/iDC1003.xml')
print(f'Original: {len(cp.reactions)} rxns, {len(cp.metabolites)} mets')

# ---- Step 1: Add transports for critical missing metabolites ----
ADD_TRANSPORTS = {
    'C00058': 'Formate',
    'C00140': 'N-Acetylglucosamine (GlcNAc)',
    'C01019': 'L-Fucose',
    'C00086': 'Urea',
}

print('\n=== Adding transport reactions ===')
added_kegg = []
existing_ids = {m.id for m in cp.metabolites}

for kegg, name in ADD_TRANSPORTS.items():
    cyto_id = f'{kegg}__cyto'
    extr_id = f'{kegg}__extr'

    if cyto_id not in existing_ids:
        print(f'  {kegg} ({name}): cyto not found, skip')
        continue
    if extr_id in existing_ids:
        print(f'  {kegg} ({name}): extr already exists, skip')
        continue

    cyto_met = cp.metabolites.get_by_id(cyto_id)
    extr_met = Metabolite(
        id=extr_id,
        formula=cyto_met.formula,
        name=cyto_met.name.replace('_C', '_extr_C'),
        compartment='C_00001',
    )
    cp.add_metabolites([extr_met])

    # Transport (reversible passive)
    tr = Reaction(id=f'T_{kegg}_cyto_extr',
                  name=f'{name} transport (cyto<->extr)',
                  lower_bound=-1000, upper_bound=1000)
    tr.add_metabolites({cyto_met: -1, extr_met: 1})
    cp.add_reactions([tr])

    # Boundary exchange
    ex = Reaction(id=f'EX_{kegg}__extr',
                  name=f'{name} exchange',
                  lower_bound=-1000, upper_bound=1000)
    ex.add_metabolites({extr_met: -1})
    cp.add_reactions([ex])

    added_kegg.append(kegg)
    print(f'  {kegg} ({name}): added ✓')

# ---- Step 2: Build CP-specific id_mapping ----
print('\n=== Building CP id_mapping ===')
ca_mapping = pd.read_csv('dat/id_mapping.csv', index_col=0)

# Direct overlaps: CP extracellular ↔ CA mapping
cp_extr_kegg = {m.id.split('__')[0]: m.id 
                for m in cp.metabolites if m.compartment == 'C_00001'}

cp_mapping_rows = []
for _, row in ca_mapping.iterrows():
    kegg = str(row['coreco'])
    if kegg in cp_extr_kegg:
        cp_mapping_rows.append({'agora': row['agora'], 'cp': cp_extr_kegg[kegg]})

# Specials: D-lactate ↔ L-lactate, NH3 ↔ NH4
if 'C00256' in cp_extr_kegg:
    cp_mapping_rows.append({'agora': 'lac_L', 'cp': cp_extr_kegg['C00256']})
    print('  Added special: C00256 (D-lactate) ↔ lac_L')
if 'C00014' in cp_extr_kegg:
    cp_mapping_rows.append({'agora': 'nh4', 'cp': cp_extr_kegg['C00014']})
    print('  Added special: C00014 (NH3) ↔ nh4')

cp_mapping = pd.DataFrame(cp_mapping_rows).drop_duplicates(subset='cp').reset_index(drop=True)
cp_mapping.to_csv('dat/id_mapping_cp.csv')
print(f'\nCP id_mapping: {len(cp_mapping)} entries → dat/id_mapping_cp.csv')

# ---- Step 3: Save modified model ----
cobra.io.write_sbml_model(cp, 'models/cp/iDC1003_processed.xml')
print(f'\nSaved: models/cp/iDC1003_processed.xml')
print(f'Final: {len(cp.reactions)} rxns, {len(cp.metabolites)} mets')

# Sanity check: model still grows?
sol = cp.optimize()
print(f'\nDefault growth: {sol.objective_value:.4f} (status: {sol.status})')

# Verify added metabolites are now mappable
print('\n=== Final coverage of OM-critical metabolites ===')
critical = {
    'Glucose': 'glc_D', 'Lactate (D→L mapping)': 'lac_L',
    'GlcNAc': 'acgam', 'Fucose L': 'fuc_L',
    'Galactose': 'gal', 'Mannose': 'man',
    'Acetate': 'ac', 'Formate': 'for',
    'Serine L': 'ser_L', 'Glycine': 'gly',
    'Aspartate L': 'asp_L', 'Glutamate L': 'glu_L',
    'NH4 (NH3 mapping)': 'nh4', 'Urea': 'urea',
    'Oxygen': 'o2',
}
mapped_agora = set(cp_mapping['agora'])
for name, agora_id in critical.items():
    status = '✓' if agora_id in mapped_agora else '✗'
    print(f'  {name:<25} ({agora_id}): {status}')
