"""AGORA metabolite ID → human-readable name mapping."""

METABOLITE_NAMES = {
    'for':   'Formate',          'succ':  'Succinate',
    'ser_L': 'L-Serine',         'fum':   'Fumarate',
    'ac':    'Acetate',          'nh4':   'Ammonium',
    'co2':   r'CO$_2$',          'h2o':   'Water',
    'o2':    r'O$_2$',           'etoh':  'Ethanol',
    'thr_L': 'L-Threonine',      'cgly':  'Cys-Gly',
    'akg':   r'α-Ketoglutarate', 'asp_L': 'L-Aspartate',
    'arg_L': 'L-Arginine',       'lac_L': 'L-Lactate',
    'mal_L': 'L-Malate',         'gln_L': 'L-Glutamine',
    'glu_L': 'L-Glutamate',      'pyr':   'Pyruvate',
    'adn':   'Adenosine',        'gsn':   'Guanosine',
    'ade':   'Adenine',          'gua':   'Guanine',
    'pi':    'Phosphate',        'amet':  'SAM',
    'amp':   'AMP',              'cit':   'Citrate',
    'gthrd': 'Glutathione',      'orn':   'L-Ornithine',
    'gly':   'Glycine',          'ala_L': 'L-Alanine',
    'his_L': 'L-Histidine',      'lys_L': 'L-Lysine',
    'met_L': 'L-Methionine',     'phe_L': 'L-Phenylalanine',
    'trp_L': 'L-Tryptophan',     'tyr_L': 'L-Tyrosine',
    'val_L': 'L-Valine',         'ile_L': 'L-Isoleucine',
    'leu_L': 'L-Leucine',        'glc_D': 'D-Glucose',
    'h2s':   r'H$_2$S',          'urea':  'Urea',
    'so4':   'Sulfate',          'hco3':  'Bicarbonate',
    'glyc':  'Glycerol',         'oaa':   'Oxaloacetate',
    'duri':  'Deoxyuridine',
    'dcyt':  'Deoxycytidine',    'dad_2': 'Deoxyadenosine',
    'dgsn':  'Deoxyguanosine',   'din':   'Deoxyinosine',
    'thymd': 'Thymidine',        'cytd':  'Cytidine',
    'uri':   'Uridine',          'ura':   'Uracil',
    'thym':  'Thymine',          'csn':   'Cytosine',
    'nac':   'Nicotinate',
}

def label(met_id):
    """Return readable name; falls back to the ID if unknown."""
    return METABOLITE_NAMES.get(met_id, met_id)
