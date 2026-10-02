"""Reference coverage and label compatibility for review, not label replacement."""
def compare(current, label, pruned, model):
    if model.startswith('Monaco') and current in {'Erythroid', 'Platelet'}:
        return 'not_covered'
    if not pruned:
        return 'uncertain'
    if model == 'HPCA_main':
        compatible = {
            'B': {'B_cell'}, 'Plasma': {'B_cell'},
            'CD4_T': {'T_cells'}, 'CD8_T': {'T_cells'},
            'NK': {'NK_cell'}, 'Monocyte_macrophage': {'Monocyte', 'Macrophage'},
            'Neutrophil': {'Neutrophils'}, 'Erythroid': {'Erythroblast'}, 'Platelet': {'Platelets'},
        }
        return 'broadly_compatible' if label in compatible[current] else 'disagreement'
    if model == 'Monaco_main':
        compatible = {
            'B': {'B cells'}, 'Plasma': {'B cells'},
            'CD4_T': {'CD4+ T cells', 'T cells'}, 'CD8_T': {'CD8+ T cells', 'T cells'},
            'NK': {'NK cells'}, 'Monocyte_macrophage': {'Monocytes'}, 'Neutrophil': {'Neutrophils'},
        }
        return 'broadly_compatible' if label in compatible[current] else 'disagreement'
    if model == 'Monaco_fine':
        if label in {'MAIT cells', 'Non-Vd2 gd T cells', 'Vd2 gd T cells'}:
            return 'T_subtype_unresolved' if current in {'CD4_T', 'CD8_T'} else 'disagreement'
        if label == 'Plasmablasts':
            return 'compatible' if current == 'Plasma' else 'disagreement'
        group = (
            'B' if label in {'Naive B cells', 'Exhausted B cells', 'Non-switched memory B cells', 'Switched memory B cells'} else
            'CD4_T' if label in {'Follicular helper T cells', 'Naive CD4 T cells', 'T regulatory cells', 'Terminal effector CD4 T cells', 'Th1 cells', 'Th1/Th17 cells', 'Th17 cells', 'Th2 cells'} else
            'CD8_T' if label in {'Central memory CD8 T cells', 'Effector memory CD8 T cells', 'Naive CD8 T cells', 'Terminal effector CD8 T cells'} else
            'NK' if label == 'Natural killer cells' else
            'Monocyte_macrophage' if label in {'Classical monocytes', 'Intermediate monocytes', 'Non classical monocytes'} else
            'Neutrophil' if label == 'Low-density neutrophils' else 'other')
        return 'compatible' if current == group else 'disagreement'
    raise ValueError(model)
