const proteins = [
  ">XP_024417592.1 transmembrane protease serine 2 [Desmodus rotundus]",
  "MALNSGLPPGVGPYYENHGFQPESLYPPQPARAPSAFSVYPAHYYPPEVPHYTVRVLTNTSTPAIRTQPK",
  "SSSGTSCTSKTKKALCITFALVAILSGAVVASVLLWKFMDSKCAGSGMECGSSGTCVSAFQWCDGILHCP",
  "SGEDENQCVRLYGPNFILQVYSPQRKSWHPVCQDNWSESYGRAACQDMGYRNSFYSSRGIADDSGATSFM",
  "KLNLSASHTDLYQQLYHSDVCSSKTVVSLRCIECGVNRKMGRQGRIVGGTSAAPGDWPWQVSLHVQNVHV",
  "CGGSIITRDWVVTAAHCVEEPLNNPRHWTAFAGILKQSSMFYGDGYRVEKVISHPNYDSKTKNNDIALMK",
  "LQTPLPFNDRVKPVCLPNPGMMLEPKQPCWISGWGATYEKGKTSDLLNAVMVPLIEPALCNHRYVYNNLI",
  "TPSMICAGYLQGNVDSCQGDSGGPLVTLKNSVWWLIGETSWGSGCAKANRPGVYGNMTVFTDWIHQQMRA",
  "NS",
  ">XP_024417600.1 transmembrane protease serine 2 [Desmodus rotundus]",
  "MALNSGLPPGVGPYYENHGFQPESLYPPQPARAPSAFSVYPAHYYPPEVPHYTVRVLTNTSTPAIRTQPK",
  "SSSGTSCTSKTKKALCITFALVAILSGAVVASVLLWKFMDSKCAGSGMECGSSGTCVSAFQWCDGILHCP",
  "SGEDENQCVRLYGPNFILQVYSPQRKSWHPVCQDNWSESYGRAACQDMGYRNSFYSSRGIADDSGATSFM",
  "KLNLSASHTDLYQQLYHSDVCSSKTVVSLRCIECGVNRKMGRQGRIVGGTSAAPGDWPWQVSLHVQNVHV",
  "CGGSIITRDWVVTAAHCVEEPLNNPRHWTAFAGILKQSSMFYGDGYRVEKVISHPNYDSKTKNNDIALMK",
  "LQTPLPFNDRVKPVCLPNPGMMLEPKQPCWISGWGATYEKGKTSDLLNAVMVPLIEPALCNHRYVYNNLI",
  "TPSMICAGYLQGNVDSCQGDSGGPLVTLKNSVWWLIGETSWGSGCAKANRPGVYGNMTVFTDWIHQQMRA",
  "NS",
  ">XP_024425698.1 angiotensin-converting enzyme 2 isoform X1 [Desmodus rotundus]",
  "MSGSSWLFLSLVAVAAAQTPTEEEARTFLENFNTEAEEWFYQNSLASWNYNTNITDENVQKMNEAEQMWS",
  "TFYERNSNIAKTYPLETIKDVNVKRQLQALQQNGLLEDKDKQLQLNAILNTMSTIYSTGKVCKPNNPQEC",
  "YLLATGLEDIMQDSKDYNERLWAWEGWRSKVGKQLRPLYEEYVVLKNEMAREKNYEDYGDYWRGDYETEG",
  "SSGYEYSRNQLIEDVENTFAEIKPLYEHLHAYVRAKLMDTYPSHISPTGCLPAHLLGDMWGRFWTNLYNL",
  "TAPFGEKPTIDVTAAMVDQSWDAQRIFKEAEKFFKSVGLFSMTQGFWDNSMLTKPDDGREVVCHPTAWDL",
  "GNKDFRIKMCTKVTMDDFLTAHHEMGHIQYDMAYANQSFLLRNGANEGFHEAVGEIMSLSVATPKHLKVL",
  "GLLPPDFHEDNETDINFLLKQALNIVGTLPFTYMLEKWRWMVFKGEIPKEQWMKKWWEMKREIVGVVEPV",
  "PHDETYCDPATLFHVANDYSFIRYYTRTIFQFQFQEALCQTAQHEGPLHKCDISNSTAAGEKLLQMLKLG",
  "KSEPWTRALENIVGKKQMDVRPLLNYFEPLFTWLKEQNRNSFVGWLTDWSPYAAESIKVRISLKSALGDK",
  "AYEWNDNEMYFFRSSIAYAMREYFSNFKNQTIPFRAEDVWVSDLKPRVSFNFFVTSPNSVSDIIPRSEVE",
  "EAIRKSRSRINDAFRLDDNSLEFLGIQPTLEPPYQPAVTIWLIAFGVVMGLVVVGIGVLIFTGIRERKRK",
  "SQETSEENPYSSMNLSKGESNPGFQNGDDVQTSF",
  ">XP_024425699.1 angiotensin-converting enzyme 2 isoform X2 [Desmodus rotundus]",
  "MWSTFYERNSNIAKTYPLETIKDVNVKRQLQALQQNGLLEDKDKQLQLNAILNTMSTIYSTGKVCKPNNP",
  "QECYLLATGLEDIMQDSKDYNERLWAWEGWRSKVGKQLRPLYEEYVVLKNEMAREKNYEDYGDYWRGDYE",
  "TEGSSGYEYSRNQLIEDVENTFAEIKPLYEHLHAYVRAKLMDTYPSHISPTGCLPAHLLGDMWGRFWTNL",
  "YNLTAPFGEKPTIDVTAAMVDQSWDAQRIFKEAEKFFKSVGLFSMTQGFWDNSMLTKPDDGREVVCHPTA",
  "WDLGNKDFRIKMCTKVTMDDFLTAHHEMGHIQYDMAYANQSFLLRNGANEGFHEAVGEIMSLSVATPKHL",
  "KVLGLLPPDFHEDNETDINFLLKQALNIVGTLPFTYMLEKWRWMVFKGEIPKEQWMKKWWEMKREIVGVV",
  "EPVPHDETYCDPATLFHVANDYSFIRYYTRTIFQFQFQEALCQTAQHEGPLHKCDISNSTAAGEKLLQML",
  "KLGKSEPWTRALENIVGKKQMDVRPLLNYFEPLFTWLKEQNRNSFVGWLTDWSPYAAESIKVRISLKSAL",
  "GDKAYEWNDNEMYFFRSSIAYAMREYFSNFKNQTIPFRAEDVWVSDLKPRVSFNFFVTSPNSVSDIIPRS",
  "EVEEAIRKSRSRINDAFRLDDNSLEFLGIQPTLEPPYQPAVTIWLIAFGVVMGLVVVGIGVLIFTGIRER",
  "KRKSQETSEENPYSSMNLSKGESNPGFQNGDDVQTSF",
  ">XP_024433200.1 LOW QUALITY PROTEIN: histo-blood group ABO system transferase [Desmodus rotundus]",
  "MVKLVRTLSQNPKCYLLYRMTFILLIALIMAFFGYCYMNPRKQRVEIIQAGASMAVRKRHDLQGVELPRM",
  "VLAQPDVLTPSRKDVLVITPWLAPIVWEGTFNIDILNEQFHLQNTTIGLTVFAIKSYIVFLKLFLETAEM",
  "YFMVGHRVNYYIFTDWPDHVPHLHLQKGRQIIILKXRSYARWEDISMHRMEMISNFFEQXFHWEVDYLVC",
  "ADVDLKFSDHVGVEILSSLFGTLHLSFFFLSRNNFPYERRPHSQAHIPEDEGDFYYTGVLFGELVMEVYR",
  "LTKTCHQVMMVDQANHIEAVWHDESHLNKYLLYHKPSKVLSPEYMWSKWLWDYSVKYQMEDLPNSVKRFR",
  "VNVLNKNQRGSN",
  ">NP_065202.2 histo-blood group ABO system transferase [Homo sapiens]",
  "MAEVLRTLAGKPKCHALRPMILFLIMLVLVLFGYGVLSPRSLMPGSLERGFCMAVREPDHLQRVSLPRMV",
  "YPQPKVLTPCRKDVLVVTPWLAPIVWEGTFNIDILNEQFRLQNTTIGLTVFAIKKYVAFLKLFLETAEKH",
  "FMVGHRVHYYVFTDQPAAVPRVTLGTGRQLSVLEVRAYKRWQDVSMRRMEMISDFCERRFLSEVDYLVCV",
  "DVDMEFRDHVGVEILTPLFGTLHPGFYGSSREAFTYERRPQSQAYIPKDEGDFYYLGGFFGGSVQEVQRL",
  "TRACHQAMMVDQANGIEAVWHDESHLNKYLLRHKPTKVLSPEYLWDQQLLGWPAVLRKLRFTAVPKNHQA",
  "VRNP",
  ">NP_001007298.1 angiotensin-converting enzyme 2 precursor [Danio rerio]",
  "MCARWLLLLALASVACCQTVEDRAREFLNKFDEEASDIMYQYTLASWAYNTDISQENADKEAEAYAIWSE",
  "YYNKMSEESNAYPIDQISDPIIKMQLQKLQDKGSGALSPDKASELRNIMSEMSTIYNTATVCKIDDPTDC",
  "QTLEPGLESIMAESRDYDERLHVWEGWRVATGMKMRPLYEKYVDLKNEAAKLNNYEDHGDYWRGDYETID",
  "DPKYSYSRDQVIEDARRIYKEILPLYKELHAYVRAKLQDVYPGHIGSDACLPAHLLGDMWGRFWTNLYPL",
  "MIPYPDRPDIDVSSAMVEQGWDEIRLFKEAEKFFMSVNMPAMFDNFWNNSMFIKPEERDVVCHPTAWDMG",
  "NRKDFRIKMCTKVNMDDFLTVHHEMGHNQYQMAYRNHPYLLRDGANEGFHEAVGEIMSLSAATPSHLQSL",
  "GLLPSDFKQDYETDINFLLKQALTIVGTLPFTYMLEEWRWQVFKAKIPKDEWMQQWWQMKRELVGVAEAV",
  "PRDETYCDPPALFHVSGDYSFIRYFTRTIYQFQFQEALCKAAGHTGPLYKCDITNSTKAGDKLRHMLELG",
  "RSMSWTRALEEVAGTTKMDSQPLLHYFSTLMEWLKEENQKNNRVPGWNVNVNPGVLTSSFINDAEISENA",
  "FKVRISLKSALGNEAYTWNANDIYLFKSTMAFAMRQYYLKEKNTDVNFTPENIHTYNETARISFKFAVMD",
  "PTKTGTVIPKAEVENAIWQERDRINGAFLLSDETLEFVGLMATLAPPKEEKITIWLVVFGVVMGVTVLAG",
  "IYLVTTGILNRKKES",
  ">XP_005169416.1 angiotensin-converting enzyme 2 isoform X2 [Danio rerio]",
  "MCARWLLLLALASVACCQTVEDRAREFLNKFDEEASDIMYQYTLASWAYNTDISQENADKEAEAYAIWSE",
  "YYNKMSEESNAYPIDQISDPIIKMQLQKLQDKGSGALSPDKASELRNIMSEMSTIYNTATVCKIDDPTDC",
  "QTLEPGLESIMAESRDYDERLHVWEGWRVATGMKMRPLYEKYVDLKNEAAKLNNYEDHGDYWRGDYETID",
  "DPKYSYSRDQVIEDARRIYKEILPLYKELHAYVRAKLQDVYPGHIGSDACLPAHLLGDMWGRFWTNLYPL",
  "MIPYPDRPDIDVSSAMVEQGWDEIRLFKEAEKFFMSVNMPAMFDNFWNNSMFIKPEERDVVCHPTAWDMG",
  "NRKDFRIKMCTKVNMDDFLTVHHEMGHNQYQMAYRNHPYLLRDGANEGFHEAVGEIMSLSAATPSHLQSL",
  "GLLPSDFKQDYETDINFLLKQALTIVGTLPFTYMLEEWRWQVFKAKIPKDEWMQQWWQMKRELVGVAEAV",
  "PRDETYCDPPALFHVSGDYSFIRYFTRTIYQFQFQEALCKAAGHTGPLYKCDITNSTKAGDKLRHMLELG",
  "RSMSWTRALEEVAGTTKMDSQPLLHYFSTLMEWLKEENQKNNRVPGWNVNVNPGVLTSSFINDAEISENA",
  "FKVRISLKSALGNEAYTWNANDIYLFKSTMAFAMRQYYLKEKNTDVNFTPENIHTYNETARISFKFAVMD",
  "PTKTGTVIPKAEVENAIWQERDRINGAFLLSDETLEFVGLMATLAPPKEEKITIWLVVFGVVMGVTVLAG",
  "IYLVTTGILNRKKKAKKAKEASVENPYDGSDDGEVNKAFEEDIEQTGL",
  ">XP_005169417.1 angiotensin-converting enzyme 2 isoform X1 [Danio rerio]",
  "MCARWLLLLALASVACCQTVEDRAREFLNKFDEEASDIMYQYTLASWAYNTDISQENADKEAEAYAIWSE",
  "YYNKMSEESNAYPIDQISDPIIKMQLQKLQDKGSGALSPDKASELRNIMSEMSTIYNTATVCKIDDPTDC",
  "QTLEPGLESIMAESRDYDERLHVWEGWRVATGMKMRPLYEKYVDLKNEAAKLNNYEDHGDYWRGDYETID",
  "DPKYSYSRDQVIEDARRIYKEILPLYKELHAYVRAKLQDVYPGHIGSDACLPAHLLGDMWGRFWTNLYPL",
  "MIPYPDRPDIDVSSAMVEQGWDEIRLFKEAEKFFMSVNMPAMFDNFWNNSMFIKPEERDVVCHPTAWDMG",
  "NRKDFRIKMCTKVNMDDFLTVHHEMGHNQYQMAYRNHPYLLRDGANEGFHEAVGEIMSLSAATPSHLQSL",
  "GLLPSDFKQDYETDINFLLKQALTIVGTLPFTYMLEEWRWQVFKAKIPKDEWMQQWWQMKRELVGVAEAV",
  "PRDETYCDPPALFHVSGDYSFIRYFTRTIYQFQFQEALCKAAGHTGPLYKCDITNSTKAGDKLRHMLELG",
  "RSMSWTRALEEVAGTTKMDSQPLLHYFSTLMEWLKEENQKNNRVPGWNVNVNPEISENAFKVRISLKSAL",
  "GNEAYTWNANDIYLFKSTMAFAMRQYYLKEKNTDVNFTPENIHTYNETARISFKFAVMDPTKTGTVIPKA",
  "EVENAIWQERDRINGAFLLSDETLEFVGLMATLAPPKEEKITIWLVVFGVVMGVTVLAGIYLVTTGILNR",
  "KKKAKKAKEASVENPYDGSDDGEVNKAFEEDIEQTGL",
  ">NP_001008623.1 transmembrane protease serine 2 [Danio rerio]",
  "MNRTETIYDNVSHVNYGFHDEEQRPPPYAPSSGVNPALSQFPAPPNTQHLPQNPHVLYPSQRYSPAHSLP",
  "HYTPQPALINTHHTVTPTIYLDTNPVNKVSRRKKWPYIVGTVVTLLVIAGVVVAVLWFYGLFACLQGRMC",
  "EADKKCVSVSLWCDGTVDCPSGEDEAQCYRLYGPNFLLQKYSNVSANWKNVCSDGFNDNLAKQACAEIGY",
  "KDTYSGYQGILSPSLDFISACQSSTAVSLKCTDCGRSTGNRIVGGTTVTSKGVWPWQVSLHYSGRHLCGG",
  "SIITPYWILTAAHCVHQFSNPGGWTVYAGYLTQSEMASASGNSVNRIVIHDFNPNTNENDIALMRLNTAL",
  "TISTNIRPVCLPNKGMSFTAQQDCYVTGWGALFSGGSSSATLQEAKIQLIDSTICNSRPVYNGLITDTMI",
  "CAGKLAGGVDSCQGDSGGPLVTNVRSLWWLLGDTSWGDGCAVRNKPGVYGNVTYFLDWIYQQMRKY",
  ">NP_056590.2 transmembrane protease serine 2 [Mus musculus]",
  "MALNSGSPPGIGPCYENHGYQSEHICPPRPPVAPNGYNLYPAQYYPSPVPQYAPRITTQASTSVIHTHPK",
  "SSGALCTSKSKKSLCLALALGTVLTGAAVAAVLLWRFWDSNCSTSEMECGSSGTCISSSLWCDGVAHCPN",
  "GEDENRCVRLYGQSFILQVYSSQRKAWYPVCQDDWSESYGRAACKDMGYKNNFYSSQGIPDQSGATSFMK",
  "LNVSSGNVDLYKKLYHSDSCSSRMVVSLRCIECGVRSVKRQSRIVGGLNASPGDWPWQVSLHVQGVHVCG",
  "GSIITPEWIVTAAHCVEEPLSSPRYWTAFAGILRQSLMFYGSRHQVEKVISHPNYDSKTKNNDIALMKLQ",
  "TPLAFNDLVKPVCLPNPGMMLDLDQECWISGWGATYEKGKTSDVLNAAMVPLIEPSKCNSKYIYNNLITP",
  "AMICAGFLQGSVDSCQGDSGGPLVTLKNGIWWLIGDTSWGSGCAKALRPGVYGNVTVFTDWIYQQMRANS",
  ">XP_006523127.1 transmembrane protease serine 2 isoform X1 [Mus musculus]",
  "MALNSGSPPGIGPCYENHGYQSEHICPPRPPVAPNGYNLYPAQYYPSPVPQYAPRITTQASTSVIHTHPK",
  "SSGALCTSKSKKSLCLALALGTVLTGAAVAAVLLWRFWDSNCSTSEMECGSSGTCISSSLWCDGVAHCPN",
  "GEDENRCVRLYGQSFILQVYSSQRKAWYPVCQDDWSESYGRAACKDMGYKNNFYSSQGIPDQSGATSFMK",
  "LNVSSGNVDLYKKLYHSDSCSSRMVVSLRCIECGVRSVKRQSRIVGGLNASPGDWPWQVSLHVQGVHVCG",
  "GSIITPEWIVTAAHCVEEPLSSPRYWTAFAGILRQSLMFYGSRHQVEKVISHPNYDSKTKNNDIALMKLQ",
  "TPLAFNDLVKPVCLPNPGMMLDLDQECWISGWGATYEKGKTSDVLNAAMVPLIEPSKCNSKYIYNNLITP",
  "AMICAGFLQGSVDSCQGDSGGPLVTLKNGIWWLIGDTSWGSGCAKALRPGVYGNVTVFTDWIYQQMRANS",
  ">XP_006523128.1 transmembrane protease serine 2 isoform X1 [Mus musculus]",
  "MALNSGSPPGIGPCYENHGYQSEHICPPRPPVAPNGYNLYPAQYYPSPVPQYAPRITTQASTSVIHTHPK",
  "SSGALCTSKSKKSLCLALALGTVLTGAAVAAVLLWRFWDSNCSTSEMECGSSGTCISSSLWCDGVAHCPN",
  "GEDENRCVRLYGQSFILQVYSSQRKAWYPVCQDDWSESYGRAACKDMGYKNNFYSSQGIPDQSGATSFMK",
  "LNVSSGNVDLYKKLYHSDSCSSRMVVSLRCIECGVRSVKRQSRIVGGLNASPGDWPWQVSLHVQGVHVCG",
  "GSIITPEWIVTAAHCVEEPLSSPRYWTAFAGILRQSLMFYGSRHQVEKVISHPNYDSKTKNNDIALMKLQ",
  "TPLAFNDLVKPVCLPNPGMMLDLDQECWISGWGATYEKGKTSDVLNAAMVPLIEPSKCNSKYIYNNLITP",
  "AMICAGFLQGSVDSCQGDSGGPLVTLKNGIWWLIGDTSWGSGCAKALRPGVYGNVTVFTDWIYQQMRANS",
  ">NP_001358344.1 angiotensin-converting enzyme 2 precursor [Homo sapiens]",
  "MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGDKWS",
  "AFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPDNPQE",
  "CLLLEPGLNEIMANSLDYNERLWAWESWRSEVGKQLRPLYEEYVVLKNEMARANHYEDYGDYWRGDYEVN",
  "GVDGYDYSRGQLIEDVEHTFEEIKPLYEHLHAYVRAKLMNAYPSYISPIGCLPAHLLGDMWGRFWTNLYS",
  "LTVPFGQKPNIDVTDAMVDQAWDAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKAVCHPTAWD",
  "LGKGDFRILMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKS",
  "IGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEP",
  "VPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNMLRL",
  "GKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYADQSIKVRISLKSALGD",
  "KAYEWNDNEMYLFRSSVAYAMRQYFLKVKNQMILFGEEDVRVANLKPRISFNFFVTAPKNVSDIIPRTEV",
  "EKAIRMSRSRINDAFRLNDNSLEFLGIQPTLGPPNQPPVSIWLIVFGVVMGVIVVGIVILIFTGIRDRKK",
  "KNKARSGENPYASIDISKGENNPGFQNTDDVQTSF",
  ">NP_068576.1 angiotensin-converting enzyme 2 precursor [Homo sapiens]",
  "MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGDKWS",
  "AFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPDNPQE",
  "CLLLEPGLNEIMANSLDYNERLWAWESWRSEVGKQLRPLYEEYVVLKNEMARANHYEDYGDYWRGDYEVN",
  "GVDGYDYSRGQLIEDVEHTFEEIKPLYEHLHAYVRAKLMNAYPSYISPIGCLPAHLLGDMWGRFWTNLYS",
  "LTVPFGQKPNIDVTDAMVDQAWDAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKAVCHPTAWD",
  "LGKGDFRILMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKS",
  "IGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEP",
  "VPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNMLRL",
  "GKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYADQSIKVRISLKSALGD",
  "KAYEWNDNEMYLFRSSVAYAMRQYFLKVKNQMILFGEEDVRVANLKPRISFNFFVTAPKNVSDIIPRTEV",
  "EKAIRMSRSRINDAFRLNDNSLEFLGIQPTLGPPNQPPVSIWLIVFGVVMGVIVVGIVILIFTGIRDRKK",
  "KNKARSGENPYASIDISKGENNPGFQNTDDVQTSF",
  ">NP_001123985.1 angiotensin-converting enzyme 2 precursor [Mus musculus]",
  "MSSSSWLLLSLVAVTTAQSLTEENAKTFLNNFNQEAEDLSYQSSLASWNYNTNITEENAQKMSEAAAKWS",
  "AFYEEQSKTAQSFSLQEIQTPIIKRQLQALQQSGSSALSADKNKQLNTILNTMSTIYSTGKVCNPKNPQE",
  "CLLLEPGLDEIMATSTDYNSRLWAWEGWRAEVGKQLRPLYEEYVVLKNEMARANNYNDYGDYWRGDYEAE",
  "GADGYNYNRNQLIEDVERTFAEIKPLYEHLHAYVRRKLMDTYPSYISPTGCLPAHLLGDMWGRFWTNLYP",
  "LTVPFAQKPNIDVTDAMMNQGWDAERIFQEAEKFFVSVGLPHMTQGFWANSMLTEPADGRKVVCHPTAWD",
  "LGHGDFRIKMCTKVTMDNFLTAHHEMGHIQYDMAYARQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKS",
  "IGLLPSDFQEDSETEINFLLKQALTIVGTLPFTYMLEKWRWMVFRGEIPKEQWMKKWWEMKREIVGVVEP",
  "LPHDETYCDPASLFHVSNDYSFIRYYTRTIYQFQFQEALCQAAKYNGSLHKCDISNSTEAGQKLLKMLSL",
  "GNSEPWTKALENVVGARNMDVKPLLNYFQPLFDWLKEQNRNSFVGWNTEWSPYADQSIKVRISLKSALGA",
  "NAYEWTNNEMFLFRSSVAYAMRKYFSIIKNQTVPFLEEDVRVSDLKPRVSFYFFVTSPQNVSDVIPRSEV",
  "EDAIRMSRGRINDVFGLNDNSLEFLGIHPTLEPPYQPPVTIWLIIFGVVMALVVVGIIILIVTGIKGRKK",
  "KNETKREENPYDSMDIGKGESNAGFQNSDDAQTSF",
  ">NP_081562.2 angiotensin-converting enzyme 2 precursor [Mus musculus]",
  "MSSSSWLLLSLVAVTTAQSLTEENAKTFLNNFNQEAEDLSYQSSLASWNYNTNITEENAQKMSEAAAKWS",
  "AFYEEQSKTAQSFSLQEIQTPIIKRQLQALQQSGSSALSADKNKQLNTILNTMSTIYSTGKVCNPKNPQE",
  "CLLLEPGLDEIMATSTDYNSRLWAWEGWRAEVGKQLRPLYEEYVVLKNEMARANNYNDYGDYWRGDYEAE",
  "GADGYNYNRNQLIEDVERTFAEIKPLYEHLHAYVRRKLMDTYPSYISPTGCLPAHLLGDMWGRFWTNLYP",
  "LTVPFAQKPNIDVTDAMMNQGWDAERIFQEAEKFFVSVGLPHMTQGFWANSMLTEPADGRKVVCHPTAWD",
  "LGHGDFRIKMCTKVTMDNFLTAHHEMGHIQYDMAYARQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKS",
  "IGLLPSDFQEDSETEINFLLKQALTIVGTLPFTYMLEKWRWMVFRGEIPKEQWMKKWWEMKREIVGVVEP",
  "LPHDETYCDPASLFHVSNDYSFIRYYTRTIYQFQFQEALCQAAKYNGSLHKCDISNSTEAGQKLLKMLSL",
  "GNSEPWTKALENVVGARNMDVKPLLNYFQPLFDWLKEQNRNSFVGWNTEWSPYADQSIKVRISLKSALGA",
  "NAYEWTNNEMFLFRSSVAYAMRKYFSIIKNQTVPFLEEDVRVSDLKPRVSFYFFVTSPQNVSDVIPRSEV",
  "EDAIRMSRGRINDVFGLNDNSLEFLGIHPTLEPPYQPPVTIWLIIFGVVMALVVVGIIILIVTGIKGRKK",
  "KNETKREENPYDSMDIGKGESNAGFQNSDDAQTSF",
  ">NP_001128571.1 transmembrane protease serine 2 isoform 1 [Homo sapiens]",
  "MPPAPPGGESGCEERGAAGHIEHSRYLSLLDAVDNSKMALNSGSPPAIGPYYENHGYQPENPYPAQPTVV",
  "PTVYEVHPAQYYPSPVPQYAPRVLTQASNPVVCTQPKSPSGTVCTSKTKKALCITLTLGTFLVGAALAAG",
  "LLWKFMGSKCSNSGIECDSSGTCINPSNWCDGVSHCPGGEDENRCVRLYGPNFILQVYSSQRKSWHPVCQ",
  "DDWNENYGRAACRDMGYKNNFYSSQGIVDDSGSTSFMKLNTSAGNVDIYKKLYHSDACSSKAVVSLRCIA",
  "CGVNLNSSRQSRIVGGESALPGAWPWQVSLHVQNVHVCGGSIITPEWIVTAAHCVEKPLNNPWHWTAFAG",
  "ILRQSFMFYGAGYQVEKVISHPNYDSKTKNNDIALMKLQKPLTFNDLVKPVCLPNPGMMLQPEQLCWISG",
  "WGATEEKGKTSEVLNAAKVLLIETQRCNSRYVYDNLITPAMICAGFLQGNVDSCQGDSGGPLVTSKNNIW",
  "WLIGDTSWGSGCAKAYRPGVYGNVMVFTDWIYRQMRADG",
  ">NP_001369649.1 transmembrane protease serine 2 isoform 3 [Homo sapiens]",
  "MALNSGSPPAIGPYYENHGYQPENPYPAQPTVVPTVYEVHPAQYYPSPVPQYAPRVLTQASNPVVCTQPK",
  "SPSGTVCTSKTKKALCITLTLGTFLVGAALAAGLLWKFMGSKCSNSGIECDSSGTCINPSNWCDGVSHCP",
  "GGEDENRCVRLYGPNFILQVYSSQRKSWHPVCQDDWNENYGRAACRDMGYKNNFYSSQGIVDDSGSTSFM",
  "KLNTSAGNVDIYKKLYHSDACSSKAVVSLRCIACGVNLNSSRQSRIVGGESALPGAWPWQVSLHVQNVHV",
  "CGGSIITPEWIVTAAHCVEKPLNNPWHWTAFAGILRQSFMFYGAGYQVEKVISHPNYDSKTKNNDIALMK",
  "LQKPLTFNDLVKPVCLPNPGMMLQPEQLCWISGWGATEEKGKTSEVLNAAKVLLIETQRCNSRYVYDNLI",
  "TPAMICAGFLQGNVDSCQGDSGGPLVTSKNNIWWLIGDTSWGSGCAKAYRPGVYGNVMVFTDWIYRQMRT",
  "ANPHGLRP",
  ">NP_005647.3 transmembrane protease serine 2 isoform 2 [Homo sapiens]",
  "MALNSGSPPAIGPYYENHGYQPENPYPAQPTVVPTVYEVHPAQYYPSPVPQYAPRVLTQASNPVVCTQPK",
  "SPSGTVCTSKTKKALCITLTLGTFLVGAALAAGLLWKFMGSKCSNSGIECDSSGTCINPSNWCDGVSHCP",
  "GGEDENRCVRLYGPNFILQVYSSQRKSWHPVCQDDWNENYGRAACRDMGYKNNFYSSQGIVDDSGSTSFM",
  "KLNTSAGNVDIYKKLYHSDACSSKAVVSLRCIACGVNLNSSRQSRIVGGESALPGAWPWQVSLHVQNVHV",
  "CGGSIITPEWIVTAAHCVEKPLNNPWHWTAFAGILRQSFMFYGAGYQVEKVISHPNYDSKTKNNDIALMK",
  "LQKPLTFNDLVKPVCLPNPGMMLQPEQLCWISGWGATEEKGKTSEVLNAAKVLLIETQRCNSRYVYDNLI",
  "TPAMICAGFLQGNVDSCQGDSGGPLVTSKNNIWWLIGDTSWGSGCAKAYRPGVYGNVMVFTDWIYRQMRA",
  "DG",
  ">NP_001129168.1 angiotensin-converting enzyme 2 precursor [Macaca mulatta]",
  "MSGSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGEKWS",
  "AFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPNNPQE",
  "CLLLDPGLNEIMEKSLDYNERLWAWEGWRSEVGKQLRPLYEEYVVLKNEMAGANHYKDYGDYWRGDYEVN",
  "GVDGYDNNRDQLIEDVERTFEEIKPLYEHLHAYVRAKLMNAYPSYISPTGCLPAHLLGDMWGRFWTNLYS",
  "LTVPFGQKPNIDVTDAMVNQAWNAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKVVCHPTAWD",
  "LGKGDFRIIMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKS",
  "IGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEP",
  "VPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLLNMLKL",
  "GESEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYADQSIKVRISLKSALGD",
  "KAYEWNDNEMYLFRSSVAYAMRTYFLEIKHQTILFGEEDVRVADLKPRISFNFYVTAPKNVSDIIPRTEV",
  "EEAIRISRSRINDAFRLNDNSLEFLGIQTTLAPPYQSPVTTWLIVFGVVMGVIVAGIVVLIFTGIRDRKK",
  "KNQARSEENPYASIDINKGENNPGFQNTDDVQTSF",
  ">XP_014982444.2 angiotensin-converting enzyme 2 isoform X1 [Macaca mulatta]",
  "MSGSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGEKWS",
  "AFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPNNPQE",
  "CLLLDPGLNEIMEKSLDYNERLWAWEGWRSEVGKQLRPLYEEYVVLKNEMARANHYKDYGDYWRGNYEVN",
  "GVDGYDYNRDQLIEDVERTFEEIKPLYEHLHAYVRAKLMNAYPSYISPTGCLPAHLLGDMWGRFWTNLYS",
  "LTVPFGQKPNIDVTDAMVNQAWNAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKVVCHPTAWD",
  "LGKGDFRIIMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKS",
  "IGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEP",
  "VPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLLNMLKL",
  "GKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYADQSIKVRISLKSALGD",
  "KAYEWNDNEMYLFRSSVAYAMRTYFLEIKHQTILFGEEDVRVADLKPRISFNFYVTAPKNVSDIIPRTEV",
  "EEAIRISRSRINDAFRLNDNSLEFLGIQTTLAPPYQSPVTTWLIVFGVVMGVIVAGIVVLIFTGIRDRKK",
  "KNQARSEENPYASIDINKGENNPGFQNTDDVQTSF",
  ">XP_028697658.1 angiotensin-converting enzyme 2 isoform X1 [Macaca mulatta]",
  "MSGSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGEKWS",
  "AFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPNNPQE",
  "CLLLDPGLNEIMEKSLDYNERLWAWEGWRSEVGKQLRPLYEEYVVLKNEMARANHYKDYGDYWRGNYEVN",
  "GVDGYDYNRDQLIEDVERTFEEIKPLYEHLHAYVRAKLMNAYPSYISPTGCLPAHLLGDMWGRFWTNLYS",
  "LTVPFGQKPNIDVTDAMVNQAWNAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKVVCHPTAWD",
  "LGKGDFRIIMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKS",
  "IGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEP",
  "VPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLLNMLKL",
  "GKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYADQSIKVRISLKSALGD",
  "KAYEWNDNEMYLFRSSVAYAMRTYFLEIKHQTILFGEEDVRVADLKPRISFNFYVTAPKNVSDIIPRTEV",
  "EEAIRISRSRINDAFRLNDNSLEFLGIQTTLAPPYQSPVTTWLIVFGVVMGVIVAGIVVLIFTGIRDRKK",
  "KNQARSEENPYASIDINKGENNPGFQNTDDVQTSF",
  ">XP_014988331.1 transmembrane protease serine 2 isoform X3 [Macaca mulatta]",
  "MALNSGSPPGVGPYYENHGYQPENPYPAQPTVAPNVYEVHPAQYYPSPVPQYTPRVLTHASNPAVCRQPK",
  "SPSGTVCTSKTKKALCVTMTLGAVLVGAALAAGLLWKFMGSKCSDSGIECDSSGTCISSSNWCDGVSHCP",
  "NGEDENRCVRLYGPNFILQVYSSQRKSWHPVCRDDWNENYARAACRDMGYKNSFYSSQGIVDNSGATSFM",
  "KLNTSAGNVDIYKKLYHSDACSSKAVVSLRCIACGVRSNLSRQSRIVGGQNALLGAWPWQVSLHVQNIHV",
  "CGGSIITPEWIVTAAHCVEKPLNSPWQWTAFVGTLRQSSMFYEKGHRVEKVISHPNYDSKTKNNDIALMK",
  "LHTPLTFNEVVKPVCLPNPGMMLEPEQHCWISGWGATQEKGKTSDVLNAAMVPLIEPRRCNNKYVYDGLI",
  "TPAMICAGFLQGTVDSCQGDSGGPLVTLKNDVWWLIGDTSWGSGCAQANRPGVYGNVTVFTDWIYRQMRA",
  "DD",
  ">XP_028701148.1 transmembrane protease serine 2 isoform X1 [Macaca mulatta]",
  "MALAGETDLDYQGEARLLLLHQAGHTEHSRYLSLLDAVDKSKMALNSGSPPGVGPYYENHGYQPENPYPA",
  "QPTVAPNVYEVHPAQYYPSPVPQYTPRVLTHASNPAVCRQPKSPSGTVCTSKTKKALCVTMTLGAVLVGA",
  "ALAAGLLWKFMGSKCSDSGIECDSSGTCISSSNWCDGVSHCPNGEDENRCVRLYGPNFILQVYSSQRKSW",
  "HPVCRDDWNENYARAACRDMGYKNSFYSSQGIVDNSGATSFMKLNTSAGNVDIYKKLYHSDACSSKAVVS",
  "LRCIACGVRSNLSRQSRIVGGQNALLGAWPWQVSLHVQNIHVCGGSIITPEWIVTAAHCVEKPLNSPWQW",
  "TAFVGTLRQSSMFYEKGHRVEKVISHPNYDSKTKNNDIALMKLHTPLTFNEVVKPVCLPNPGMMLEPEQH",
  "CWISGWGATQEKGKTSDVLNAAMVPLIEPRRCNNKYVYDGLITPAMICAGFLQGTVDSCQGDSGGPLVTL",
  "KNDVWWLIGDTSWGSGCAQANRPGVYGNVTVFTDWIYRQMRADD",
  ">XP_028701149.1 transmembrane protease serine 2 isoform X2 [Macaca mulatta]",
  "MSGLRILSLHKADIGHTEHSRYLSLLDAVDKSKMALNSGSPPGVGPYYENHGYQPENPYPAQPTVAPNVY",
  "EVHPAQYYPSPVPQYTPRVLTHASNPAVCRQPKSPSGTVCTSKTKKALCVTMTLGAVLVGAALAAGLLWK",
  "FMGSKCSDSGIECDSSGTCISSSNWCDGVSHCPNGEDENRCVRLYGPNFILQVYSSQRKSWHPVCRDDWN",
  "ENYARAACRDMGYKNSFYSSQGIVDNSGATSFMKLNTSAGNVDIYKKLYHSDACSSKAVVSLRCIACGVR",
  "SNLSRQSRIVGGQNALLGAWPWQVSLHVQNIHVCGGSIITPEWIVTAAHCVEKPLNSPWQWTAFVGTLRQ",
  "SSMFYEKGHRVEKVISHPNYDSKTKNNDIALMKLHTPLTFNEVVKPVCLPNPGMMLEPEQHCWISGWGAT",
  "QEKGKTSDVLNAAMVPLIEPRRCNNKYVYDGLITPAMICAGFLQGTVDSCQGDSGGPLVTLKNDVWWLIG",
  "DTSWGSGCAQANRPGVYGNVTVFTDWIYRQMRADD",
  ">XP_014971938.2 histo-blood group ABO system transferase isoform X1 [Macaca mulatta]",
  "MGARARRRRPFLAGVPGDPRPPPAPRCPHPSPASASPGGGRRAARTGRGRLELTPPQGCRAEGGSGDQTL",
  "SHGGAAADADCTAFREPDGLQRVSLPRMVYPQPKVLTPCRKDVLVVTPWLAPIVWEGTFNIDILNEQFRL",
  "QNTTIGLTVFAIKKYVAFLKLFLETAEKHFMVGHRVHYYVFTDQPAAVPRVALGTGRQLSVLGVRAYKRW",
  "QDVSMRRMEMISDFCERRFLSEVDYLVCADVDMEFRDHVGVEILTPLFGTLHPAFYGSSREAFTYERRPQ",
  "SQAYIPKDEGDFYYMGAFFGGSVQEVQRLTRACHQAMMVDQANSIEAVWHDESHLNKYLLRHKPTKVLSP",
  "EYLWDQQLLGWPAVLRKLRFVAVPKNHQAVRNP",
  ">XP_014971939.2 histo-blood group ABO system transferase isoform X3 [Macaca mulatta]",
  "MGARARRRRPFLAGVPGDPRPPPAPRCPHPSPASASPGGGRRAARTGRGRLELTPPQGCRAEGGSGDQTL",
  "SHGGAAADADWMVYPQPKVLTPCRKDVLVVTPWLAPIVWEGTFNIDILNEQFRLQNTTIGLTVFAIKKYV",
  "AFLKLFLETAEKHFMVGHRVHYYVFTDQPAAVPRVALGTGRQLSVLGVRAYKRWQDVSMRRMEMISDFCE",
  "RRFLSEVDYLVCADVDMEFRDHVGVEILTPLFGTLHPAFYGSSREAFTYERRPQSQAYIPKDEGDFYYMG",
  "AFFGGSVQEVQRLTRACHQAMMVDQANSIEAVWHDESHLNKYLLRHKPTKVLSPEYLWDQQLLGWPAVLR",
  "KLRFVAVPKNHQAVRNP",
  ">XP_014971940.2 histo-blood group ABO system transferase isoform X4 [Macaca mulatta]",
  "MAELLQTLTGKPKCYTFLPMILFLMMLVLVLFGYGVLSPRSLMPGSLERGFCTAFREPDGLQRVSLPRMV",
  "YPQPKVLTPCRKDVLVVTPWLAPIVWEGTFNIDILNEQFRLQNTTIGLTVFAIKKYVAFLKLFLETAEKH",
  "FMVGHRVHYYVFTDQPAAVPRVALGTGRQLSVLGVRAYKRWQDVSMRRMEMISDFCERRFLSEVDYLVCA",
  "DVDMEFRDHVGVEILTPLFGTLHPAFYGSSREAFTYERRPQSQAYIPKDEGDFYYMGAFFGGSVQEVQRL",
  "TRACHQAMMVDQANSIEAVWHDESHLNKYLLRHKPTKVLSPEYLWDQQLLGWPAVLRKLRFVAVPKNHQA",
  "VRNP",
  ">XP_014971941.1 histo-blood group ABO system transferase isoform X6 [Macaca mulatta]",
  "MAELLQTLTGKPKCYTFLPMILFLMMLVLVLFGYGVLSPRSLMPGSLERGFWMVYPQPKVLTPCRKDVLV",
  "VTPWLAPIVWEGTFNIDILNEQFRLQNTTIGLTVFAIKKYVAFLKLFLETAEKHFMVGHRVHYYVFTDQP",
  "AAVPRVALGTGRQLSVLGVRAYKRWQDVSMRRMEMISDFCERRFLSEVDYLVCADVDMEFRDHVGVEILT",
  "PLFGTLHPAFYGSSREAFTYERRPQSQAYIPKDEGDFYYMGAFFGGSVQEVQRLTRACHQAMMVDQANSI",
  "EAVWHDESHLNKYLLRHKPTKVLSPEYLWDQQLLGWPAVLRKLRFVAVPKNHQAVRNP",
  ">XP_028690481.1 histo-blood group ABO system transferase isoform X2 [Macaca mulatta]",
  "MGARARRRRPFLAGVPGDPRPPPAPRCPHPSPASASPGGGRRAARTGRGRLELTPPQGCRAEGGSGDQTL",
  "SHGGAAADADCTAFREPDGLQRVSLPRMVYPQPKVLTPWKDVLVVTPWLAPIVWEGTFNIDILNEQFRLQ",
  "NTTIGLTVFAIKKYVAFLKLFLETAEKHFMVGHRVHYYVFTDQPAAVPRVALGTGRQLSVLGVRAYKRWQ",
  "DVSMRRMEMISDFCERRFLSEVDYLVCADVDMEFRDHVGVEILTPLFGTLHPAFYGSSREAFTYERRPQS",
  "QAYIPKDEGDFYYMGAFFGGSVQEVQRLTRACHQAMMVDQANSIEAVWHDESHLNKYLLRHKPTKVLSPE",
  "YLWDQQLLGWPAVLRKLRFVAVPKNHQAVRNP",
  ">XP_028690482.1 histo-blood group ABO system transferase isoform X5 [Macaca mulatta]",
  "MAELLQTLTGKPKCYTFLPMILFLMMLVLVLFGYGVLSPRSLMPGSLERGFCTAFREPDGLQRVSLPRMV",
  "YPQPKVLTPWKDVLVVTPWLAPIVWEGTFNIDILNEQFRLQNTTIGLTVFAIKKYVAFLKLFLETAEKHF",
  "MVGHRVHYYVFTDQPAAVPRVALGTGRQLSVLGVRAYKRWQDVSMRRMEMISDFCERRFLSEVDYLVCAD",
  "VDMEFRDHVGVEILTPLFGTLHPAFYGSSREAFTYERRPQSQAYIPKDEGDFYYMGAFFGGSVQEVQRLT",
  "RACHQAMMVDQANSIEAVWHDESHLNKYLLRHKPTKVLSPEYLWDQQLLGWPAVLRKLRFVAVPKNHQAV",
  "RNP",
  ">XP_028690483.1 histo-blood group ABO system transferase isoform X7 [Macaca mulatta]",
  "MAELLQTLTGKPKCYTFLPMILFLMMLVLVLFGYGVLSPRSLMPGSLERGFWMVYPQPKVLTPWKDVLVV",
  "TPWLAPIVWEGTFNIDILNEQFRLQNTTIGLTVFAIKKYVAFLKLFLETAEKHFMVGHRVHYYVFTDQPA",
  "AVPRVALGTGRQLSVLGVRAYKRWQDVSMRRMEMISDFCERRFLSEVDYLVCADVDMEFRDHVGVEILTP",
  "LFGTLHPAFYGSSREAFTYERRPQSQAYIPKDEGDFYYMGAFFGGSVQEVQRLTRACHQAMMVDQANSIE",
  "AVWHDESHLNKYLLRHKPTKVLSPEYLWDQQLLGWPAVLRKLRFVAVPKNHQAVRNP",
  ">XP_028690484.1 histo-blood group ABO system transferase isoform X8 [Macaca mulatta]",
  "MVYPQPKVLTPCRKDVLVVTPWLAPIVWEGTFNIDILNEQFRLQNTTIGLTVFAIKKYVAFLKLFLETAE",
  "KHFMVGHRVHYYVFTDQPAAVPRVALGTGRQLSVLGVRAYKRWQDVSMRRMEMISDFCERRFLSEVDYLV",
  "CADVDMEFRDHVGVEILTPLFGTLHPAFYGSSREAFTYERRPQSQAYIPKDEGDFYYMGAFFGGSVQEVQ",
  "RLTRACHQAMMVDQANSIEAVWHDESHLNKYLLRHKPTKVLSPEYLWDQQLLGWPAVLRKLRFVAVPKNH",
  "QAVRNP",
  ">XP_028690485.1 histo-blood group ABO system transferase isoform X9 [Macaca mulatta]",
  "MVYPQPKVLTPWKDVLVVTPWLAPIVWEGTFNIDILNEQFRLQNTTIGLTVFAIKKYVAFLKLFLETAEK",
  "HFMVGHRVHYYVFTDQPAAVPRVALGTGRQLSVLGVRAYKRWQDVSMRRMEMISDFCERRFLSEVDYLVC",
  "ADVDMEFRDHVGVEILTPLFGTLHPAFYGSSREAFTYERRPQSQAYIPKDEGDFYYMGAFFGGSVQEVQR",
  "LTRACHQAMMVDQANSIEAVWHDESHLNKYLLRHKPTKVLSPEYLWDQQLLGWPAVLRKLRFVAVPKNHQ",
  "AVRNP",
  ">NP_001277373.1 histo-blood group ABO system transferase isoform b [Mus musculus]",
  "MNLRGRPKCNFLHLGILPFAVFVLVFFGYLFLSFRSQNLGHPGAVTRNAYLQPRVLKPTYVVFLKLFLET",
  "AEQHFMVGHKVIYYVFTDRPADVPQVILGAGRQLVVLTVRNYTRWQDVSMHRMEMISHFSERRFLREVDY",
  "LVCADADMKFSDHVGVEILSTFFGTLHPGFYSSSREAFTYERRPQSQAYIPWDRGDFYYGGAFFGGSVLE",
  "VYHLTKACHEAMMEDKANGIEPVWHDESYLNKYLLYHKPTKVLSPEYLWDQQLLGWPSIMKKLRYVAVPK",
  "DHQAIRN",
  ">NP_109643.3 histo-blood group ABO system transferase isoform a [Mus musculus]",
  "MNLRGRPKCNFLHLGILPFAVFVLVFFGYLFLSFRSQNLGHPGAVTRNAYLQPRVLKPTRKDVLVLTPWL",
  "APIIWEGTFNIDILNEQFRIRNTTIGLTVFAIKKYVVFLKLFLETAEQHFMVGHKVIYYVFTDRPADVPQ",
  "VILGAGRQLVVLTVRNYTRWQDVSMHRMEMISHFSERRFLREVDYLVCADADMKFSDHVGVEILSTFFGT",
  "LHPGFYSSSREAFTYERRPQSQAYIPWDRGDFYYGGAFFGGSVLEVYHLTKACHEAMMEDKANGIEPVWH",
  "DESYLNKYLLYHKPTKVLSPEYLWDQQLLGWPSIMKKLRYVAVPKDHQAIRN",
  ">XP_006498519.1 histo-blood group ABO system transferase isoform X1 [Mus musculus]",
  "MRQEPEQVRQSDDRSAVGPIRMGAQSHASSWKSFWFQQGTGSYTQNENRGREKDRGDPCRSGRPKCNFLH",
  "LGILPFAVFVLVFFGYLFLSFRSQNLGHPGAVTRNAYLQPRVLKPTRKDVLVLTPWLAPIIWEGTFNIDI",
  "LNEQFRIRNTTIGLTVFAIKKYVVFLKLFLETAEQHFMVGHKVIYYVFTDRPADVPQVILGAGRQLVVLT",
  "VRNYTRWQDVSMHRMEMISHFSERRFLREVDYLVCADADMKFSDHVGVEILSTFFGTLHPGFYSSSREAF",
  "TYERRPQSQAYIPWDRGDFYYGGAFFGGSVLEVYHLTKACHEAMMEDKANGIEPVWHDESYLNKYLLYHK",
  "PTKVLSPEYLWDQQLLGWPSIMKKLRYVAVPKDHQAIRN",
  ">XP_006498521.1 histo-blood group ABO system transferase isoform X2 [Mus musculus]",
  "MRQEPEQVRQSDDRSAVGYVVFLKLFLETAEQHFMVGHKVIYYVFTDRPADVPQVILGAGRQLVVLTVRN",
  "YTRWQDVSMHRMEMISHFSERRFLREVDYLVCADADMKFSDHVGVEILSTFFGTLHPGFYSSSREAFTYE",
  "RRPQSQAYIPWDRGDFYYGGAFFGGSVLEVYHLTKACHEAMMEDKANGIEPVWHDESYLNKYLLYHKPTK",
  "VLSPEYLWDQQLLGWPSIMKKLRYVAVPKDHQAIRN",
];
