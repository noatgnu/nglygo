package configuration

type Configuration struct {
	RootPath string
	CPUCore int
	SpeciesFile string
	OutputFolder string
	UniprotTabulatedFile string
	QueryFastaFile string
	OutputBlastFile string
	BlastDB string
	DefaultTree string
	PythonPath string
	FastTreePath string
	ClustalOmegaPath string
	BlastDBCMDPath string
	MakeBlastDBPath string
	BlastpPath string
	CodeMLPath string
	PhyMLPath string
	RaxMLPath string
	ProcessTreePythonScriptPath string
	PruneTreePythonScripthPath string
	MakeBlastDBInputFasta string
	MakeBlastDBOutputFilteredFasta string
	ServerPort int
	BlastTargetNumber int
	BlastEvalue float64
	BlastIdentityThreshold float64
	BlastCoverageThreshold float64
	PhylogeneticTreeConstructionSpeed string
	ASRProgramPath string
	ASRMatrixPath string
	HomologLengthMaxProportion int
}