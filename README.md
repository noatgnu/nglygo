# ancestral
 
## Configuration


For operation of the compiled library, 

```
-c    Create sequon DB
-config string
Configuration file path (default "config.json")
-m    Create Blastp DB
-w    Activate webserver
```

These parameters can be used with a `configuration.json` or `config.json` file in order create a BlastP database or N-linked sequon database.

```json
{
  "CPUCore": 7,
  "SpeciesFile": "C:\\Users\\localadmin\\PycharmProjects\\ancestralplay\\spike_protein_blast.uniprot.txt",
  "OutputFolder": "D:\\GoProject\\ancestral\\result\\spike.protein.20210428",
  "UniprotTabulatedFile": "C:\\Users\\localadmin\\Downloads\\spike.protein.blast.uniprot.tab",
  "QueryFastaFile": "C:\\Users\\localadmin\\PycharmProjects\\ancestralplay\\spike.fasta",
  "OutputBlastFile": "C:\\Users\\localadmin\\PycharmProjects\\ancestralplay\\spike_protein.20210428.blast.tsv",
  "BlastDB": "D:\\GoProject\\ancestral\\spike_protein_blast.uniprot_customDB",
  "BlastTargetNumber": 1200,
  "BlastEValue": 0.001,
  "BlastIdentityThreshold": 0,
  "BlastCoverageThreshold": 0,
  "PhylogeneticTreeConstructionSpeed": "fast",
  "DefaultTree": "",
  "RootPath": "D:\\GoProject\\ancestral\\result",
  "PythonPath": "C:\\Program Files\\Anaconda3\\python.exe",
  "FastTreePath": "D:\\GoProject\\ancestral\\FastTree.exe",
  "PhyMLPath": "D:\\PhyML-3.1\\PhyML-3.1_win32.exe",
  "ClustalOmegaPath": "clustalo",
  "BlastDBCMDPath": "blastdbcmd",
  "MakeBlastDBPath": "makeblastdb",
  "BlastpPath": "blastp",
  "CodeMLPath": "D:\\paml4.9e\\bin\\codeml.exe",
  "RaxMLPath": "raxML",
  "ProcessTreePythonScriptPath": "D:\\GoProject\\ancestral\\src\\github.com\\noatgnu\\ancestral\\processtree.py",
  "PruneTreePythonScripthPath": "D:\\GoProject\\ancestral\\src\\github.com\\noatgnu\\ancestral\\prunetree.py",
  "ASRProgramPath": "D:\\paml4.9e\\bin\\codeml.exe",
  "ASRMatrixPath" : "D:\\paml4.9e\\dat\\wag.dat",

  "MakeBlastDBInputFasta": "C:\\Users\\localadmin\\PycharmProjects\\ancestralplay\\spike_protein_blast.uniprot.fasta",
  "MakeBlastDBOutputFilteredFasta": "C:\\Users\\localadmin\\PycharmProjects\\ancestralplay\\spike_protein_blast.uniprot.fasta_filtered",
  "ServerPort": 8080

}
```

For `-m` parameter to work, `MakeBlastDBInputFasta`, `MakeBlastDBOutputFilteredFasta`, `MakeBlastDBPath`,  and `BlastDB` have to be defined in the configuration json.

For `-s` parameter to work:

1. `CPUCore` >= 1
2. `SpeciesFile` is a text find where the name of the species could be found separated by newlines.
3. `OutputFolder` output folder for the analysis. Each subfolder one level down from this folder is the result from each query of the Blastp operation.
4. `UniprotTabulatedFile` tabulated file download from Uniprot containing Uniprot accession ID, species name, entry name, sequence, and topological domain

