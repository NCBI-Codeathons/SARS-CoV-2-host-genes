#!/usr/bin/env python3

TAXID_TO_ABBR = {
    # Manual added abbreviations.
    9430:'Vba',
    # Abbreviations from NCBI internal database: taxpropdump gpipe_abbr
    2711:'Csi', 2850:'Ptc', 2903:'Ehu', 3055:'Cre', 3067:'Vca', 3218:'Ppa', 
    3330:'Pgl', 3332:'Psi', 3352:'Pta', 3369:'Cjn', 3635:'Ghi', 3641:'Tcc', 
    3649:'Cpp', 3659:'Cst', 3694:'Pth', 3695:'Ptd', 3702:'At', 3708:'Bna', 
    3711:'Bra', 3712:'Bol', 3726:'Rsa', 3750:'Mdo', 3760:'Ppe', 3818:'Ahy', 
    3827:'Car', 3847:'Gma', 3880:'Mtr', 3885:'Pvu', 3886:'Pco', 3917:'Vun', 
    3983:'Mes', 4072:'Can', 4081:'Les', 4097:'Nta', 4111:'Sml', 4113:'Stu', 
    4155:'Mgu', 4232:'Han', 4236:'Lsa', 4513:'Hv', 4530:'Os', 4533:'Obr', 
    4547:'Sof', 4558:'Sbi', 4565:'Ta', 4577:'Zm', 4608:'Fpr', 4679:'Ace', 
    4754:'Pca', 4787:'Pin', 4896:'Spo', 4911:'Kma', 4914:'Kwa', 4927:'Wan', 
    4929:'Pgu', 4931:'Sba', 4932:'Sce', 4934:'Skl', 4952:'Yli', 4956:'Zro', 
    4959:'Dha', 5037:'Aca', 5039:'Ade', 5057:'Acv', 5061:'Ani', 5067:'Apa', 
    5141:'Ncr', 5207:'Fne', 5270:'Uma', 5272:'Mvi', 5306:'Pch', 5334:'Sco', 
    5346:'Cci', 5476:'Cal', 5478:'Cgl', 5482:'Ctr', 5501:'Cim', 5518:'Gze', 
    5664:'Lma', 5671:'Lin', 5691:'Tbr', 5692:'Tco', 5693:'Tcr', 5699:'Tvi', 
    5722:'Tva', 5741:'Gin', 5759:'Ehi', 5791:'Ppo', 5807:'Cpa', 5811:'Tgo', 
    5821:'Pbe', 5825:'Pcb', 5833:'Pfa', 5849:'Pga', 5850:'Pkn', 5854:'Pre', 
    5855:'Pvi', 5861:'Pyo', 5865:'Bbo', 5874:'Tan', 5875:'Tpa', 5888:'Pte', 
    5911:'Tth', 6035:'Ecu', 6087:'Hma', 6182:'Sja', 6183:'Sma', 6238:'Cbr', 
    6239:'Cel', 6279:'Bma', 6305:'Mha', 6334:'Tsp', 6376:'Apo', 6412:'Hro', 
    6500:'Acl', 6526:'Bgl', 6584:'Sso', 6613:'Esc', 6669:'Dpu', 6689:'Lva', 
    6945:'Isc', 7029:'Aps', 7070:'Tca', 7091:'Bmo', 7102:'Hvi', 7159:'Aae', 
    7165:'Aga', 7175:'Cps', 7176:'Cpi', 7213:'Ccp', 7217:'Dan', 7220:'Der', 
    7222:'Dgr', 7227:'Dm', 7230:'Dmo', 7234:'Dpe', 7237:'Dps', 7238:'Dse', 
    7240:'Dsi', 7244:'Dvi', 7245:'Dya', 7260:'Dwi', 7370:'Mds', 7394:'Gms', 
    7425:'Nvi', 7426:'Ngi', 7427:'Nlo', 7460:'Ame', 7462:'Ado', 7463:'Afl', 
    7656:'Pli', 7668:'Spu', 7719:'Cin', 7739:'Bfl', 7757:'Pma', 7868:'Cmi', 
    7897:'Lch', 7918:'Loc', 7955:'Dr', 7998:'Ipu', 8018:'Oke', 8019:'Oki', 
    8022:'Omy', 8023:'One', 8030:'Ssa', 8049:'Gmr', 8078:'Fhe', 8083:'Xma', 
    8090:'Ola', 8128:'Oni', 8153:'Hbu', 8355:'Xl', 8364:'Str', 8478:'Cpt', 
    8479:'Cpc', 8496:'Ami', 8839:'Apl', 8954:'Fpe', 9031:'Gga', 9103:'Mga', 
    9258:'Oan', 9302:'Smc', 9305:'Sha', 9315:'Meu', 9337:'Tvu', 9361:'Dno', 
    9365:'Eeu', 9368:'Aal', 9371:'Ete', 9417:'Aja', 9447:'Lca', 9483:'Cja', 
    9509:'Age', 9519:'Lla', 9523:'Cmo', 9534:'Cae', 9541:'Mfa', 9544:'Mmu', 
    9545:'Mne', 9555:'Pan', 9556:'Pcy', 9587:'Hkl', 9593:'Ggo', 9597:'Ppn', 
    9598:'Ptr', 9600:'Ppy', 9601:'Pab', 9606:'Hs', 9614:'Cla', 9615:'Cfa', 
    9646:'Aml', 9668:'Mpu', 9669:'Mpt', 9685:'Fca', 9708:'Oro', 9733:'Oor', 
    9739:'Ttn', 9785:'Laf', 9796:'Eca', 9823:'Ssc', 9901:'Bbi', 9913:'Bt', 
    9925:'Chi', 9940:'Oar', 9978:'Opr', 9986:'Ocu', 10029:'Cgr', 10036:'Mau', 
    10042:'Pmn', 10090:'Mm', 10096:'Msp', 10116:'Rn', 10141:'Cpr', 10160:'Ode', 
    10181:'Hgb', 10224:'Sko', 10228:'Tad', 12957:'Acp', 13146:'Mun', 13616:'Mdm', 
    13735:'Psn', 15368:'Bdi', 27288:'Nca', 27291:'Spa', 27293:'Sse', 27679:'Sbo', 
    27828:'Con', 28377:'Acr', 28737:'Eed', 28985:'Kla', 29170:'Acn', 29730:'Gra', 
    29760:'Vvi', 29850:'Ggr', 29883:'Lbi', 30195:'Bte', 30286:'Mte', 30301:'Bsc', 
    30538:'Vpa', 30608:'Mmr', 30611:'Oga', 31033:'Tru', 31136:'Pbi', 31234:'Crm', 
    31870:'Ggm', 33085:'Ein', 33169:'Ego', 33178:'Ate', 33548:'Cgu', 34305:'Lja', 
    34317:'Egl', 34358:'Sex', 34638:'Moc', 34765:'Odi', 34839:'Cln', 35128:'Tps', 
    35608:'Aan', 35883:'Ini', 36630:'Nfi', 36914:'Lel', 37293:'Ana', 37690:'Ptf', 
    38033:'Cgb', 38654:'Asi', 38727:'Pvr', 39638:'Aph', 39758:'Mde', 40233:'Cpe', 
    40559:'Bfu', 42254:'Sar', 42374:'Cdu', 43179:'Itr', 44015:'Gto', 44689:'Ddi', 
    45157:'Cme', 45351:'Nve', 45683:'Tla', 46681:'Edi', 46689:'Fhy', 47247:'Lco', 
    47664:'Ptp', 49390:'Cca', 51029:'Hgl', 51337:'Jja', 51453:'Hje', 51511:'Csa', 
    51620:'Kya', 52253:'Cso', 53485:'Pts', 54126:'Ppf', 55529:'Gth', 57918:'Fve', 
    58336:'Svo', 59463:'Mlu', 59479:'Rfe', 59538:'Pho', 59689:'Aly', 59729:'Tgu', 
    59894:'Fal', 60710:'Cpy', 60711:'Csb', 61452:'Nne', 61853:'Nle', 62324:'Afn', 
    64495:'Ror', 66913:'Ifu', 67593:'Pso', 69293:'Gac', 70448:'Ota', 71647:'Ppi', 
    72004:'Bgr', 72036:'Lsl', 73337:'Csm', 73824:'Pba', 74940:'Ots', 77166:'Dpo', 
    78454:'Sla', 79327:'Sme', 79684:'Moh', 81526:'Mov', 85681:'Ccl', 88036:'Smo', 
    88211:'Pci', 89462:'Bbu', 90988:'Ppr', 94289:'Shi', 99883:'Tni', 104421:'Cfl', 
    106582:'Mze', 106881:'Lme', 109478:'Mbr', 109996:'Rra', 114524:'Sku', 
    114525:'Smi', 117187:'Gmo', 120017:'Uho', 121225:'Phu', 127582:'Tma', 
    128810:'Shc', 130081:'Gsu', 132113:'Bim', 135651:'Csp', 143302:'Ccr', 
    143995:'Mro', 148305:'Mgr', 162425:'Eni', 164328:'Pra', 169999:'Pme', 
    170000:'Ppc', 199306:'Cpo', 225164:'Lgi', 227086:'Bnt', 237895:'Cho', 
    242159:'Osp', 246437:'Tch', 278021:'Alo', 281687:'Cjp', 303518:'Pny', 
    338618:'Afp', 381046:'Kth', 400682:'Aqu', 419612:'Cfe', 464988:'Had', 
    529818:'Ttr', 554065:'Cva', 586133:'Npa', 610380:'Hsa', 746128:'Afu', 
    2052682:'Pul', 2587412:'Pas'
}


def get_short_taxname(taxid, name=None):
    if taxid in TAXID_TO_ABBR:
        return TAXID_TO_ABBR[taxid]
    if name:
        return name.split()[-1]
    return 'unknown'


def get_short_taxname_from_gene(gene_data):
    name = gene_data.get('taxname', None)
    name = gene_data.get('commonName', None)
    tax_id = int(gene_data['taxId'])
    return get_short_taxname(tax_id, name)
