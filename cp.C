void cp(){
	TGrid::Connect("alien:///");
	TFile::Cp("JDL2", "alien:///alice/cern.ch/user/k/kimb/jetjt/JDL2");
	//TFile::Cp("pythiaChargedDijet.C", "alien:///alice/cern.ch/user/k/kimb/jetjt/pythiaChargedDijet.C");
	//TFile::Cp("pythia8230.tar", "alien:///alice/cern.ch/user/k/kimb/jetjt/pythia8230.tar");
	TFile::Cp("Makefile", "alien:///alice/cern.ch/user/k/kimb/jetjt/Makefile");
	//TFile::Cp("src.tar", "alien:///alice/cern.ch/user/k/kimb/jetjt/src.tar");
	//TFile::Cp("pythia_config_Monash2013.cmnd", "alien:///alice/cern.ch/user/k/kimb/jetjt/pythia_config_Monash2013.cmnd");
	//TFile::Cp("cardAlice_pp.input", "alien:///alice/cern.ch/user/k/kimb/jetjt/cardAlice_pp.input");
	//TFile::Cp("pythia_config_4C.cmnd", "alien:///alice/cern.ch/user/k/kimb/jetjt/pythia_config_4C.cmnd");
	TFile::Cp("fastpythia", "alien:///alice/cern.ch/user/k/kimb/bin/fastpythia");
}
