"""Project-wide matplotlib font/style settings.
Importing this module applies Arial-equivalent sans-serif font and
embeds fonts properly in PDF/PS outputs (Type 42 / TrueType)."""

import matplotlib as mpl

# Font: Arial first, then open-source fallbacks (Linux/WSL 환경 호환)
mpl.rcParams['font.family']     = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']

# PDF/PS에 폰트 임베드 (저널 submission requirement 대응)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42

# Math text도 같은 sans-serif 계열로 통일 (μ, ×, subscript 등)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm']      = 'Arial'
mpl.rcParams['mathtext.it']      = 'Arial:italic'
mpl.rcParams['mathtext.bf']      = 'Arial:bold'
