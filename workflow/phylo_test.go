package workflow

import (
	"testing"
	"encoding/json"
	"github.com/noatgnu/ancestral/alignio"
	"os"
	"encoding/csv"
)

var blastMapJson =  `{"id":{"1":"\u003eNP_001004690.1 olfactory receptor 2M5 [Homo sapiens] \u003elcl|A3KFT3.1 RecName: Full=Olfactory receptor 2M5 \u003elcl|EAW77216.1 olfactory receptor, family 2, subfamily M, member 5 [Homo sapiens]","10":"\u003eXP_013846920.1 olfactory receptor 2M3 [Sus scrofa]","11":"\u003eXP_008983880.2 PREDICTED: olfactory receptor 2M3-like [Callithrix jacchus] \u003elcl|XP_017822200.1 PREDICTED: olfactory receptor 2M3-like [Callithrix jacchus]","12":"\u003eXP_013974368.1 PREDICTED: olfactory receptor 2M3-like [Canis lupus familiaris]","13":"\u003eXP_011281372.1 PREDICTED: olfactory receptor 2M3 [Felis catus]","14":"\u003eXP_001493933.2 PREDICTED: olfactory receptor 2M3-like [Equus caballus]","15":"\u003eXP_004008798.1 PREDICTED: olfactory receptor 2M3-like [Ovis aries] \u003elcl|XP_012024540.1 PREDICTED: olfactory receptor 2M3-like [Ovis aries musimon]","16":"\u003eXP_013820894.1 PREDICTED: olfactory receptor 2M3-like [Capra hircus]","17":"\u003eXP_002689149.1 PREDICTED: olfactory receptor 2M3 isoform X2 [Bos taurus] \u003elcl|XP_588102.2 PREDICTED: olfactory receptor 2M3 isoform X2 [Bos taurus] \u003elcl|DAA27557.1 TPA: olfactory receptor, family 2, subfamily B, member 2-like [Bos taurus]","18":"\u003eXP_003267291.2 PREDICTED: olfactory receptor 2M4 [Nomascus leucogenys]","19":"\u003eXP_005368980.1 PREDICTED: olfactory receptor 2M3-like [Microtus ochrogaster]","2":"\u003eXP_016807162.1 PREDICTED: olfactory receptor 2M5 [Pan troglodytes]","20":"\u003eEDL77975.1 olfactory receptor 1570 (predicted) [Rattus norvegicus]","21":"\u003eXP_017194966.1 PREDICTED: olfactory receptor 2M3-like [Oryctolagus cuniculus]","22":"\u003eNP_666662.1 olfactory receptor 164 [Mus musculus] \u003elcl|XP_021040176.1 olfactory receptor 2M3-like [Mus caroli] \u003elcl|AAL61310.1 olfactory receptor MOR279-2 [Mus musculus] \u003elcl|EDK97548.1 mCG1037977 [Mus musculus]","23":"\u003eXP_021065526.1 olfactory receptor 2M3-like [Mus pahari]","24":"\u003eXP_019840928.1 PREDICTED: olfactory receptor 2M2-like [Bos indicus]","25":"\u003eXP_021032412.1 olfactory receptor 56-like [Mus caroli]","26":"\u003eXP_001372591.2 PREDICTED: olfactory receptor 56-like [Monodelphis domestica]","27":"\u003eXP_002809238.1 PREDICTED: olfactory receptor 2T12-like [Pongo abelii]","28":"\u003eXP_001511658.2 PREDICTED: olfactory receptor 2T33-like [Ornithorhynchus anatinus]","3":"\u003eXP_018894989.1 PREDICTED: olfactory receptor 2M7 [Gorilla gorilla gorilla]","4":"\u003eXP_003815363.1 PREDICTED: olfactory receptor 2M7 [Pan paniscus]","5":"\u003eXP_003893682.1 PREDICTED: olfactory receptor 2M3 [Papio anubis]","6":"\u003eXP_001094325.2 PREDICTED: olfactory receptor 2M3 [Macaca mulatta]","7":"\u003eEHH61754.1 hypothetical protein EGM_19841 [Macaca fascicularis]","8":"\u003eXP_007988300.1 PREDICTED: olfactory receptor 2M3-like [Chlorocebus sabaeus]","9":"\u003eXP_012613559.1 olfactory receptor 2M3 [Microcebus murinus]"},"organism":{"1":"Homo sapiens","10":"Sus scrofa","11":"Callithrix jacchus","12":"Canis lupus familiaris","13":"Felis catus","14":"Equus caballus","15":"Ovis aries","16":"Capra hircus","17":"Bos taurus","18":"Nomascus leucogenys","19":"Microtus ochrogaster","2":"Pan troglodytes","20":"Rattus norvegicus","21":"Oryctolagus cuniculus","22":"Mus musculus","23":"Mus pahari","24":"Bos indicus","25":"Mus caroli","26":"Monodelphis domestica","27":"Pongo abelii","28":"Ornithorhynchus anatinus","3":"Gorilla gorilla gorilla","4":"Pan paniscus","5":"Papio anubis","6":"Macaca mulatta","7":"Macaca fascicularis","8":"Chlorocebus sabaeus","9":"Microcebus murinus"},"accession":{"1":"\u003eNP_001004690.1","10":"\u003eXP_013846920.1","11":"\u003eXP_008983880.2","12":"\u003eXP_013974368.1","13":"\u003eXP_011281372.1","14":"\u003eXP_001493933.2","15":"\u003eXP_004008798.1","16":"\u003eXP_013820894.1","17":"\u003eXP_002689149.1","18":"\u003eXP_003267291.2","19":"\u003eXP_005368980.1","2":"\u003eXP_016807162.1","20":"\u003eEDL77975.1","21":"\u003eXP_017194966.1","22":"\u003eNP_666662.1","23":"\u003eXP_021065526.1","24":"\u003eXP_019840928.1","25":"\u003eXP_021032412.1","26":"\u003eXP_001372591.2","27":"\u003eXP_002809238.1","28":"\u003eXP_001511658.2","3":"\u003eXP_018894989.1","4":"\u003eXP_003815363.1","5":"\u003eXP_003893682.1","6":"\u003eXP_001094325.2","7":"\u003eEHH61754.1","8":"\u003eXP_007988300.1","9":"\u003eXP_012613559.1"},"filename":"\"D:\\\\GoProject\\\\ancestral\\\\result\\\\test2\\\\A3KFT3\\\\compiled.filtered.fasta\"","MatchSourceID":"1","SourceSeq":{"Id":"A3KFT3","Seq":"MAWENQTFNSDFILLGIFNHSPTHTFLFFLVLAIFSVAFMGNSVMVLLIYLDTQLHTPMYFLLSQLFLMDLMLICSTVPKMAFNYLSGSKSISMAGCATQIFFYVSLLGSECFLLAVMSYDRYIAICHPLRYTNLMRPKICGLMTAFSWILGSMDAIIDAVATFSFSYCGSREIAHFFCDFPSLLILSCNDTSIFEKVLFICCIVMIVFPVAIIIASYARVILAVIHMGSGEGRRKAFTTCSSHLMVVGMYYGAGLFMYIRPTSDRSPMQDKLVSVFYTILTPMLNPLIYSLRNKEVTRALRKVLGKGKCGE","Length":312,"Species":"Homo sapiens (Human)","Name":"OR2M5_HUMAN","TopDomain":[{"Start":1,"Stop":25,"Type":"Extracellular"},{"Start":50,"Stop":57,"Type":"Cytoplasmic"},{"Start":80,"Stop":100,"Type":"Extracellular"},{"Start":121,"Stop":139,"Type":"Cytoplasmic"},{"Start":159,"Stop":195,"Type":"Extracellular"},{"Start":220,"Stop":236,"Type":"Cytoplasmic"},{"Start":260,"Stop":272,"Type":"Extracellular"},{"Start":293,"Stop":311,"Type":"Cytoplasmic"}]}}
`

func TestRemapTopDom(t *testing.T) {
	b := BlastMap{}
	err := json.Unmarshal([]byte(blastMapJson), &b)
	if err != nil {
		t.Error(err)
	}
	alignment := alignio.ReadPhylip(`compiled.phy`)
	RemapTopDom(alignment.Alignment[b.MatchSourceID], b.SourceSeq.TopDomain, "test.txt")
	f, err := os.Open("test.txt")
	if err != nil {
		t.Error(err)
	}
	reader := csv.NewReader(f)
	reader.Comma = '\t'
	rows, err := reader.ReadAll()
	if err != nil {
		t.Error(err)
	}
	l := len(rows)
	if len(rows) != 9 {
		t.Errorf("Expected %v, got %v", 9, l)
	}
	f.Close()
	err = os.Remove("test.txt")
	if err != nil {
		t.Error(err)
	}
}
