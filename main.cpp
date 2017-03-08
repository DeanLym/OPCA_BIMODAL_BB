/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
#include "opca_bimodal.h"
#include "sim_model.hpp"
#include "generate_opca_model.hpp"
#include "util_funs.hpp"


int main ( int argc , char ** argv ) {

	// display:
	int dim_ = 80;
	if (argc < 2){
		cout << "Please specify input file." << std::endl;
	}
	ifstream ifs;
	ifs.open(argv[1]);
	vector<double> xi;
	double temp;
	for(int i=0;i<dim_;i++){
		ifs >> temp;
		xi.push_back(temp);
	}
	TranformUniform2Normal(dim_,xi);
	OPCA_BIMODAL* opca_bm_  = GenerateOPCAModel();
	double Nc = 3600;
	vector<double> m;
	m.resize(Nc);
	bool opca_flag;
	opca_flag = opca_bm_->GenerateOPCARealization(&(xi[0]) , &(m[0]));
	if(!opca_flag) // If O-PCA returns false, the solution is wrong
	{
		string log_file("opca_err_log");
		LogOpcaError(log_file.c_str(), xi, dim_);
		return true;       // the evaluation failed
	}
	vector<double> perm;
	perm.resize(Nc);
	GeneratePerm(Nc, &(m[0]), &(perm[0]));

#ifdef DEBUG
	SaveData("xi.debug",dim_,&(xi[0]));
	SaveData("m.debug",Nc,&(m[0]));
	SaveData("perm.debug",Nc,&(perm[0]));
#endif

	SimCtrl* sim = GetSimulationModel(&perm[0]);

	sim->display_level_ = 0;
	//cout << "Begin to run simulation." <<endl;
	sim->RunSim();
	//cout << "Finished running simulation." << endl;
#if defined(DEBUG) || defined(PRED)
	sim->OutputResult();
#endif
	sim->hm_ = new CHistoryMatching;
	ifstream in;
	in.open("HIST_FILE.DATA");
	string hist_file, temp_str;
	if(in.is_open())
		in >> temp_str >> hist_file;
	else{
		throw runtime_error("Can not open HIST_FILE.DATA");
		//			cout << "Can not open HIST_FILE.DATA" << endl;
	}
	in.close();
	//		cout << hist_file << endl;
	sim->hm_->SetHMTarget(hist_file.c_str());
	//cout << "Finished setting history matching target." << endl;
	//		cout << "Set history matching target..." << endl;
	vector<double> d = sim->hm_->GetData(sim->std_well_);

	SaveData("d.txt",d.size(),&(d[0]));
	delete opca_bm_;
	return true;       // the evaluation succeeded
}
