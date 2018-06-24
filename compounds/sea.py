__author__ = 'baiqing'
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
import os
import math
import cPickle
from glaucoma import settings

TARGET_FOLDER_BASE = os.path.join(settings.BASE_DIR, 'target-pkl')
TARGET_LIST = ['CHEMBL1827',
               'CHEMBL205',
               'CHEMBL2034',
               'CHEMBL210',
               'CHEMBL211',
               'CHEMBL213',
               'CHEMBL214',
               'CHEMBL1867',
               'CHEMBL216',
               'CHEMBL217',
               'CHEMBL1951',
               'CHEMBL2056',
               'CHEMBL220',
               'CHEMBL221',
               'CHEMBL222',
               'CHEMBL223',
               'CHEMBL224',
               'CHEMBL225',
               'CHEMBL226',
               'CHEMBL229',
               'CHEMBL230',
               'CHEMBL232',
               'CHEMBL1987',
               'CHEMBL1881',
               'CHEMBL2015',
               'CHEMBL1940',
               'CHEMBL241',
               'CHEMBL1783',
               'CHEMBL1821',
               'CHEMBL2035',
               'CHEMBL1942',
               'CHEMBL1916',
               'CHEMBL245',
               'CHEMBL1904',
               'CHEMBL246',
               'CHEMBL1833',
               'CHEMBL1900',
               'CHEMBL1785',
               'CHEMBL252',
               'CHEMBL251',
               'CHEMBL1836',
               'CHEMBL254',
               'CHEMBL256',
               'CHEMBL4716',
               'CHEMBL2885',
               'CHEMBL3421',
               'CHEMBL3025',
               'CHEMBL4884',
               'CHEMBL261',
               'CHEMBL3710',
               'CHEMBL1914',
               'CHEMBL4425',
               'CHEMBL2652',
               'CHEMBL275',
               'CHEMBL3231',
               'CHEMBL3805',
               'CHEMBL2609',
               'CHEMBL3969',
               'CHEMBL4789',
               'CHEMBL2326',
               'CHEMBL3510',
               'CHEMBL2973',
               'CHEMBL4619',
               'CHEMBL286',
               'CHEMBL288',
               'CHEMBL3977',
               'CHEMBL290',
               'CHEMBL291',
               'CHEMBL1980',
               'CHEMBL4235',
               'CHEMBL3878',
               'CHEMBL3119',
               'CHEMBL3012',
               'CHEMBL3242',
               'CHEMBL2717',
               'CHEMBL3729',
               'CHEMBL3912',
               'CHEMBL3594',
               'CHEMBL2072',
               'CHEMBL4296',
               'CHEMBL4187',
               'CHEMBL5163',
               'CHEMBL5451',
               'CHEMBL3535',
               'CHEMBL3746',
               'CHEMBL3836',
               'CHEMBL4640',
               'CHEMBL4408',
               'CHEMBL4409',
               'CHEMBL4267',
               'CHEMBL5932']
FP_PARAM = {
    'topological_hashed': {
        "mean": 0.000026936641120031898,
        "sd": 0.00398734520007183,
        "sd_exp": 0.5145149536573593,
        "tc": 0.66,
        "fp_func": lambda m: AllChem.GetHashedTopologicalTorsionFingerprint(m)
    },
    'atompair_hashed': {
        "mean": 0.00002541280101763038,
        "sd": 0.003864080986890242,
        "sd_exp": 0.5118760617288712,
        "tc": 0.61,
        "fp_func": lambda m: AllChem.GetHashedAtomPairFingerprint(m)
    },
    'maccs': {
        "mean": 0.00003686608558399457,
        "sd": 0.00484043909504593,
        "sd_exp": 0.5199608769322853,
        "tc": 0.88,
        "fp_func": lambda m: AllChem.GetMACCSKeysFingerprint(m)
    },
    'morgan_hashed': {
        "mean": 0.00002318732496154415,
        "sd": 0.003829515743309366,
        "sd_exp": 0.5101931872278405,
        "tc": 0.65,
        "fp_func": lambda m: AllChem.GetHashedMorganFingerprint(m, 2)
    }
}


def raw_score(target_mol_pkl, mol_fp, cutoff):
    # try:
    sim_list = list()
    for el_fp in target_mol_pkl:
        sim = TanimotoSimilarity(mol_fp, el_fp)
        if sim >= cutoff:
            sim_list.append(sim)
    rawscore = sum(sim_list)
    if rawscore > 0:
        return round(rawscore, 3)
    return None


def z_score(rs, size, mean, sd, sd_exp):
    return (rs - size * mean) / (sd * size ** sd_exp)


def p_value(z):
    x = -math.exp(-z * math.pi / math.sqrt(6) - 0.577215665)
    if z > 28:
        return -(x + x ** 2 / 2 + x ** 3 / 6)
    else:
        return 1 - math.exp(x)


def pred2_build_in_function(smiles, target_list=TARGET_LIST):
    result = list()
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        # todo: invalid rdkit molecule
        print 'invalid smiles'
        return None
    # for fp_name, fp_parm in FP_PARAM.iteritems():
    #     mol_fp = fp_parm['fp_func'](mol)
    #     print fp_name
    #     fp_result = dict()
    fp_dict = dict()
    for fp_name, fp_param in FP_PARAM.iteritems():
        fp_dict[fp_name] = fp_param['fp_func'](mol)

    for idx, chembl_id in enumerate(target_list):
        # print idx
        target_result = dict()
        target_result['smiles'] = smiles
        for fp_name, fp_param in FP_PARAM.iteritems():
            mol_fp = fp_dict[fp_name]
            target_mol_pkl = cPickle.load(open(os.path.join(TARGET_FOLDER_BASE, fp_name, chembl_id), 'r'))
            rs = raw_score(target_mol_pkl, mol_fp, fp_param['tc'])
            if rs:
                zscore = z_score(rs, len(target_mol_pkl), fp_param['mean'], fp_param['sd'], fp_param['sd_exp'])
                pvalue = p_value(zscore)
                target_result[fp_name] = pvalue

        if target_result:
            target_result['chembl_id'] = chembl_id
            target_result['smile'] = smiles
            result.append(target_result)
    return result

def pred2_all_target_list(smiles, target_list=TARGET_LIST):
    result = list()
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        # todo: invalid rdkit molecule
        print 'invalid smiles'
        return None
    # for fp_name, fp_parm in FP_PARAM.iteritems():
    #     mol_fp = fp_parm['fp_func'](mol)
    #     print fp_name
    #     fp_result = dict()
    fp_dict = dict()
    for fp_name, fp_param in FP_PARAM.iteritems():
        fp_dict[fp_name] = fp_param['fp_func'](mol)

    for idx, chembl_id in enumerate(target_list):
        # print idx
        target_result = dict()
        target_result['smiles'] = smiles
        target_result['chembl_id'] = chembl_id
        for fp_name, fp_param in FP_PARAM.iteritems():
            mol_fp = fp_dict[fp_name]
            target_mol_pkl = cPickle.load(open(os.path.join(TARGET_FOLDER_BASE, fp_name, chembl_id), 'r'))
            rs = raw_score(target_mol_pkl, mol_fp, fp_param['tc'])
            if rs:
                zscore = z_score(rs, len(target_mol_pkl), fp_param['mean'], fp_param['sd'], fp_param['sd_exp'])
                pvalue = p_value(zscore)
                target_result[fp_name] = pvalue
            else:
                target_result[fp_name] = None
        result.append(target_result)

        # if target_result:
        #     target_result['chembl_id'] = chembl_id
        #     # target_result['smile'] = smiles
        #     result.append(target_result)
    return result

def pred2_single_target(smiles, target):
    result = list()
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        # todo: invalid rdkit molecule
        print 'invalid smiles'
        return None
    # for fp_name, fp_parm in FP_PARAM.iteritems():
    #     mol_fp = fp_parm['fp_func'](mol)
    #     print fp_name
    #     fp_result = dict()
    fp_dict = dict()
    for fp_name, fp_param in FP_PARAM.iteritems():
        fp_dict[fp_name] = fp_param['fp_func'](mol)

    for idx, chembl_id in enumerate([target]):
        # print idx
        target_result = dict()
        target_result['smiles'] = smiles
        target_result['chembl_id'] = chembl_id
        for fp_name, fp_param in FP_PARAM.iteritems():
            mol_fp = fp_dict[fp_name]
            target_mol_pkl = cPickle.load(open(os.path.join(TARGET_FOLDER_BASE, fp_name, chembl_id), 'r'))
            rs = raw_score(target_mol_pkl, mol_fp, fp_param['tc'])
            if rs:
                zscore = z_score(rs, len(target_mol_pkl), fp_param['mean'], fp_param['sd'], fp_param['sd_exp'])
                pvalue = p_value(zscore)
                target_result[fp_name] = pvalue
            else:
                target_result[fp_name] = None
        result.append(target_result)
        # if target_result:
        #     target_result['chembl_id'] = chembl_id
        #     target_result['smile'] = smiles
        #     result.append(target_result)
    return result


def pred(smiles, target_list):
    result = dict()
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        # todo: invalid rdkit molecule
        print 'invalid smiles'
        return None
    for fp_name, fp_parm in FP_PARAM.iteritems():
        mol_fp = fp_parm['fp_func'](mol)
        print fp_name
        fp_result = dict()
        for idx, chembl_id in enumerate(target_list):
            # print idx
            target_mol_pkl = cPickle.load(open(os.path.join(TARGET_FOLDER_BASE, fp_name, chembl_id), 'r'))
            rs = raw_score(target_mol_pkl, mol_fp, fp_parm['tc'])
            if rs:
                zscore = z_score(rs, len(target_mol_pkl), fp_parm['mean'], fp_parm['sd'], fp_parm['sd_exp'])
                pvalue = p_value(zscore)
                fp_result[chembl_id] = {
                    'p-value': pvalue
                }
        result[fp_name] = fp_result
    return result


def _test():
    # target_list = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))[:100]
    # TARGET_LIST = ['CHEMBL2034',
    #                'CHEMBL4267',
    #                'CHEMBL2717',
    #                'CHEMBL1987',
    #                'CHEMBL5932',
    #                'CHEMBL286',
    #                'CHEMBL3119']
    smiles = 'CC1CCCN(C1C)C(=O)c2csc(Nc3ccc(C)cc3)n2'
    print pred2(smiles, TARGET_LIST)


# if __name__ == '__main__':
#     _test()
