"""CP-specific lumen integration helpers."""
import cobra
from cobra import Metabolite, Reaction

def cp_join_lumen_rxns(cp_model, mapping_df):
    """Add lumen exchanges to CP model (CP extracellular ↔ shared lumen).

    Args:
        cp_model: processed CP cobrapy Model (extracellular = C_00001)
        mapping_df: DataFrame with columns 'agora' and 'cp'

    Returns:
        Modified cp_model with lumen ('u') compartment + exchange reactions.
    """
    lumen_mets = {}
    existing_rxn_ids = {r.id for r in cp_model.reactions}
    rxns_to_add = []
    mets_to_add = []

    for _, row in mapping_df.iterrows():
        agora = row['agora']
        cp_id = row['cp']

        try:
            cp_met = cp_model.metabolites.get_by_id(cp_id)
        except KeyError:
            continue

        # Lumen metabolite
        lumen_id = f'{agora}_u'
        if agora not in lumen_mets:
            try:
                lumen_met = cp_model.metabolites.get_by_id(lumen_id)
            except KeyError:
                lumen_met = Metabolite(id=lumen_id, compartment='u',
                                       name=f'{agora} (lumen)',
                                       formula=cp_met.formula or '')
                mets_to_add.append(lumen_met)
            lumen_mets[agora] = lumen_met
        else:
            lumen_met = lumen_mets[agora]

        # CP extr ↔ lumen exchange
        new_rxn_id = f'Ex_{cp_id}'
        if new_rxn_id not in existing_rxn_ids:
            r = Reaction(id=new_rxn_id,
                         name=f'{cp_id} <-> lumen',
                         lower_bound=-1000, upper_bound=1000)
            r.add_metabolites({cp_met: -1, lumen_met: 1})
            rxns_to_add.append(r)
            existing_rxn_ids.add(new_rxn_id)

        # Disable original CP boundary (EX_<cp_id>__extr)
        original_ex_id = f'EX_{cp_id}'
        try:
            orig = cp_model.reactions.get_by_id(original_ex_id)
            if orig.id != new_rxn_id and len(orig.metabolites) == 1:
                orig.lower_bound = 0
                orig.upper_bound = 0
        except KeyError:
            pass

    # System boundary: EX_u_<agora>
    for agora, lumen_met in lumen_mets.items():
        sys_id = f'EX_u_{agora}'
        if sys_id not in existing_rxn_ids:
            r = Reaction(id=sys_id, name=f'EX {agora} (lumen->env)',
                         lower_bound=-1000, upper_bound=1000)
            r.add_metabolites({lumen_met: -1})
            rxns_to_add.append(r)

    cp_model.add_metabolites(mets_to_add)
    cp_model.add_reactions(rxns_to_add)
    return cp_model


def categorize_joint_rxn(reaction):
    """Categorize a reaction in joint model by compartment signature.
    Handles both CA ('cytosol') and CP ('C_00001') as fungal compartments.
    """
    comps = set(m.compartment for m in reaction.metabolites)
    fungal_comps = {'cytosol', 'C_00001'}
    if comps & fungal_comps and 'u' in comps and 'e' not in comps:
        return 'fungal'
    if comps == {'e', 'u'}: return 'BAC'
    if comps == {'u'}: return 'system'
    return None


def make_joint_model_cp(cp_model, bac_model, second_lumen_ext_rxns,
                         second_lumen_ext_mets, diet, copy_module, source):
    """CP-aware joint model. Mirhakkak make_joint_model와 차이점:
    1. Fungal biomass: 'e_Biomass__cyto' (CA의 'scerbiomasspseudoreaction' 대신)
    2. O2 bound: diet 따름 (강제 anaerobic 제거 — OM context)
    """
    joint_model = copy_module.deepcopy(cp_model)
    joint_model.solver = 'glpk'
    for i in range(len(second_lumen_ext_rxns)):
        joint_model.add_reactions([second_lumen_ext_rxns[i]])
        second_lumen_ext_rxns[i].add_metabolites({second_lumen_ext_mets[i]: -1})
        second_lumen_ext_rxns[i].bounds = (0, 1000)
        if second_lumen_ext_rxns[i].id in diet.keys():
            second_lumen_ext_rxns[i].lower_bound = diet[second_lumen_ext_rxns[i].id]
    joint_model.add_reactions(bac_model.reactions)

    # O2 from diet (NOT hardcoded anaerobic)
    if 'EX_u_o2' in diet.keys():
        joint_model.reactions.get_by_id('EX_u_o2').lower_bound = diet['EX_u_o2']

    bac_biomass_id = 'model2_biomass' if source == 'agora' else 'model2_Growth'
    biomass_id_second_model = next(
        r.id for r in joint_model.reactions if r.id.startswith(bac_biomass_id))

    biomass1 = joint_model.reactions.get_by_id('e_Biomass__cyto')
    biomass2 = joint_model.reactions.get_by_id(biomass_id_second_model)
    biomass1.lower_bound = 0.001
    biomass2.lower_bound = 0.001
    joint_model.objective = [biomass1, biomass2]
    return joint_model
