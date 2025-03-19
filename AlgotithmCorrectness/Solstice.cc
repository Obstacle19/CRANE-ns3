#include <iostream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <unistd.h> 

using namespace std;

int used_ports;


struct OSwitchDemandMatrixElem {
	int row;
	int col;
	uint64_t value;
};

class OSwitchDemandMatrix {
public:
	int Numrow;
	int Numcol;
	uint64_t** m_array;

	OSwitchDemandMatrix(int row, int col){
		Numrow = row;
		Numcol = col;
		m_array = new uint64_t* [Numrow];
		for (int i = 0;i < Numrow; i++){
			m_array[i] = new uint64_t[Numcol];
		}
		for (int i = 0; i < Numrow; i++){
			for (int j = 0; j < Numcol; j++){
				m_array[i][j] = 0;
			}
		}
	}
	~OSwitchDemandMatrix(){
		for (int i = 0;i < Numrow; i++){
			delete[] m_array[i];
		}
		delete[] m_array;
	}

	OSwitchDemandMatrixElem* Find_Max(){
		uint64_t _max_value = 0;
		OSwitchDemandMatrixElem* __res = new OSwitchDemandMatrixElem;
		__res->col = __res->row = 0;
		__res->value = 0;
		for (int i = 0; i < Numrow; i++)
		{
			for (int j = 0; j < Numcol; j++)
			{
				if (m_array[i][j] > _max_value){
					__res->row = i;
					__res->col = j;
					__res->value = m_array[i][j];
					_max_value = m_array[i][j];
				}
			}
		}
		return __res;
	}

	void Set_matrix(uint64_t** mtx, int rnum, int cnum){
		for(int i = 0;i < rnum; i++){
			for (int j = 0;j < cnum; j++){
				m_array[i][j] = mtx[i][j];
			}
		}
	}
};

//-----------------------------------Hungary Algorithm---------------------------------------------------------
bool bpm(int u, std::vector<std::vector<int>>& bpGraph, std::vector<bool>& seen, std::vector<int>& matchR) {
    for (int v = 0; v < bpGraph[0].size(); v++) {
        if (bpGraph[u][v] && !seen[v]) {
            seen[v] = true;
            if (matchR[v] < 0 || bpm(matchR[v], bpGraph, seen, matchR)) {
                matchR[v] = u;
                return true;
            }
        }
    }
    return false;
}

int maxBPM(std::vector<std::vector<int>>& bpGraph, std::vector<int>& matchR) {
    int U = bpGraph.size();
    int V = bpGraph[0].size();
    for (int u = 0; u < U; u++) {
        std::vector<bool> seen(V, false);
        bpm(u, bpGraph, seen, matchR);
    }
    
    int maxMatches = 0;
    for (int v = 0; v < V; v++) {
        if (matchR[v] != -1) {
            maxMatches++;
        }
    }
    return maxMatches;
}

int** MyOswtichHungrayAlgorithm(int** cmatrix, int nodenums){
	int U = nodenums;
	int V = nodenums;
	std::vector<std::vector<int>> bpGraph(U, std::vector<int>(V, 0));
	std::vector<int> matchR(V, -1);
    for (int i = 0; i < U; i++) {
        for (int j = 0; j < V; j++) {
            bpGraph[i][j] = cmatrix[i][j];
        }
    }
	int maxMatch = maxBPM(bpGraph, matchR);
	int** result;
	result = new int*[nodenums];
	for (int i = 0; i < nodenums; i++){
		result[i] = new int[nodenums];
	}
	for (int i = 0; i < nodenums; i++){
		for (int j = 0;j < nodenums; j ++){
			result[i][j] = 0;
		}
	}
	for (int v = 0; v < V; v++) {
        if (matchR[v] != -1) {
            result[matchR[v]][v] = 1;
        }
    }
	return result;
}

//-----------------------------------Hungary Algorithm End---------------------------------------------------------

class Solstice_output{	//The return value of Solstice Algorithm
public:
	int Numcol;
	int Numrow;
	int duration_nums; //The 'm', means there are m circuit configurations.
	std::vector<uint64_t> duration_time;	//Store each configure's duration time.
	bool*** circuit_configures;	//circuit_configures[durationnum][sender][reciever]
	uint64_t** PSwitchMatrix;


	void Init(){
		PSwitchMatrix = new uint64_t* [Numrow];
		for (int i = 0;i < Numrow; i++){
			PSwitchMatrix[i] = new uint64_t[Numcol];
		}
	}

	void cleanup(){
		for (int i = 0;i < duration_nums; i++){
			for (int j = 0; j < Numrow; j++){
				delete[] circuit_configures[i][j];
			}
			delete[] circuit_configures[i];
		}
		delete[] circuit_configures;
		for (int i = 0;i < Numrow; i++){
			delete[] PSwitchMatrix[i];
		}
		delete[] PSwitchMatrix;
	}

	void Set_duration_time (int* durations, int n){
		for (int i = 0;i < n; i++){
			duration_time[i] = durations[i];
		}
	}
	void Set_circuit_configures(bool value, int dura_time, int rnum, int cnum){
		circuit_configures[dura_time][rnum][cnum] = value;
	}
	void Set_PSwitchMatrix(uint64_t** PSM, int rnum, int cnum){
		for(int i = 0;i < rnum; i++){
			for (int j = 0;j < cnum; j++){
				PSwitchMatrix[i][j] = PSM[i][j];
			}
		}
	}


	Solstice_output(int Numcol, int Numrow, int duration_nums){
		this->duration_nums = duration_nums;
		this->Numcol = Numcol;
		this->Numrow = Numrow;
		Init();
	}

	~Solstice_output(){
		cleanup();
	}

	void display() {
		cout << "There are " << duration_nums << " duration nums:\n";
		for (int i = 0; i < duration_nums; i++){
			cout << "The number " << i << " st/nd/th duation is:\n";
			for (int a = 0; a < Numrow; a++){
				for(int b = 0; b < Numcol; b++){
					cout<< circuit_configures[i][a][b] << ' ';
				}
				cout << '\n';
			}
			cout << "This configure's duration is " << duration_time[i] << " ns.\n";
		}
		
		cout << "And the Packet demand matrix is:\n";
		for (int a = 0; a < Numrow; a++){
			for(int b = 0; b < Numcol; b++){
				cout<< PSwitchMatrix[a][b] << ' ';
			}
			cout << '\n';
		}
	}

};


class SliceOutput {
public:
	bool valid;	//A succesful slice output or not.
	int rowNum;
	int colNum;
	bool** Permutation_Matrix;

	SliceOutput(bool va, int rnum, int cnum){
		valid = va;
		rowNum = rnum;
		colNum = cnum;
		Permutation_Matrix = new bool* [rnum];
		for (int i = 0;i < rnum; i++){
			Permutation_Matrix[i] = new bool[cnum];
		}
		for (int i = 0;i < rnum; i++){
			for (int j = 0;j < cnum; j++){
				Permutation_Matrix[i][j] = 0;
			}
		}
	}
	~SliceOutput(){
		for (int i = 0;i < rowNum; i++){
			delete[] Permutation_Matrix[i];
		}
		delete[] Permutation_Matrix;
	}

};

SliceOutput* BigSlice(OSwitchDemandMatrix* matrix, uint64_t threshold){
	SliceOutput* __res;
	//Return value, default to be not valid, determine it will be valid or not later.
	__res = new SliceOutput(false, used_ports, used_ports);

	//prepare the matrix
	int** tempmatrix = new int*[used_ports];
	for (int i = 0;i < used_ports; i++){
		tempmatrix[i] = new int[used_ports];
	}
	// for (int i = 0;i < used_ports; i++){
	// 	for (int j = 0;j < used_ports; j++){
	// 		tempmatrix[i][j] = matrix->m_array[i][j];
	// 	}
	// }
//-------------------------------------------------------------------------
	
	//First Do D'' <- ZeroEntiresBelow(D', r) AND B-<BinaryMatrixOf(D'')
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			if (matrix->m_array[i][j] < threshold){
				tempmatrix[i][j] = 0;
			}
			else{
				tempmatrix[i][j] = 1;
			}
		}
	}
	//Now find the perfect matching of B, we used hungary algorithm to find it.
	int** matchresult = MyOswtichHungrayAlgorithm(tempmatrix, used_ports);
	int resultsum = 0;
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			resultsum += matchresult[i][j];
		}
	}
	if (resultsum == used_ports){
		__res->valid = 1;
		for (int i = 0; i < used_ports; i++){
			for (int j = 0; j < used_ports; j++){
				__res->Permutation_Matrix[i][j] = matchresult[i][j] > 0;
			}
		}
	}

	for (int i = 0;i < used_ports; i++){
		delete[] tempmatrix[i];
	}
	delete[] tempmatrix;
	return __res;
}

OSwitchDemandMatrix* QuickStuff(OSwitchDemandMatrix* matrix){
	uint64_t Rsums[used_ports];//Ri
	uint64_t Csums[used_ports];//Ci
	uint64_t diameter;//Fi = max of(Ri, Ci)
	OSwitchDemandMatrix* __res;
	__res = new OSwitchDemandMatrix(used_ports, used_ports);
	uint64_t** tempmatrix = new uint64_t*[used_ports];
	for (int i = 0;i < used_ports; i++){
		tempmatrix[i] = new uint64_t[used_ports];
	}
	for (int i = 0;i < used_ports; i++){
		for (int j = 0;j < used_ports; j++){
			tempmatrix[i][j] = matrix->m_array[i][j];
		}
	}
	__res->Set_matrix(tempmatrix, used_ports, used_ports);
	for (int i = 0;i < used_ports; i++){
		delete[] tempmatrix[i];
	}
	delete[] tempmatrix;
//-----------------------------------------------------------------------------------------
	//caculate Ri and Ci
	for (int i = 0;i < used_ports; i++){
		Rsums[i] = 0;
		for (int j = 0; j < used_ports; j++)
		{
			Rsums[i] += __res->m_array[i][j];
		}
	}
	for (int j = 0;j < used_ports; j++){
		Csums[j] = 0;
		for (int i = 0;i < used_ports; i++){
			Csums[j] += __res->m_array[i][j];
		}
	}
	// Find diameter
	diameter = 0;
	for (int i = 0; i < used_ports; i++){
		if(Rsums[i] > diameter){
			diameter = Rsums[i];
		}
		if(Csums[i] > diameter){
			diameter = Csums[i];
		}
	}
	
	uint64_t temp;
	//Ajust the non-zero values
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			if (__res->m_array[i][j] > 0) {
				temp = diameter - std::max(Rsums[i], Csums[j]);
				__res->m_array[i][j] += temp;
				Rsums[i] += temp;
				Csums[j] += temp;
			}
		}
	}
	//Ajust the zero values
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			if (__res->m_array[i][j] == 0) {
				temp = diameter - std::max(Rsums[i], Csums[j]);
				__res->m_array[i][j] += temp;
				Rsums[i] += temp;
				Csums[j] += temp;
			}
		}
	}
	return __res;
}

void SolsticeUpdateDiameter(OSwitchDemandMatrix* DemandM, uint64_t* Diameter){
	uint64_t Rsums[used_ports];//Ri
	uint64_t Csums[used_ports];//Ci
	*Diameter = 0;
	for (int i = 0;i < used_ports; i++){
		Rsums[i] = 0;
		for (int j = 0; j < used_ports; j++)
		{
			Rsums[i] += DemandM->m_array[i][j];
		}
	}
	for (int j = 0;j < used_ports; j++){
		Csums[j] = 0;
		for (int i = 0;i < used_ports; i++){
			Csums[j] += DemandM->m_array[i][j];
		}
	}
	for (int i = 0; i < used_ports; i++){
		if(Rsums[i] > *Diameter){
			cout << "Diameter changed, before:" << *Diameter;
			*Diameter = Rsums[i];
			cout << " after " << *Diameter << '\n';
		}
		if(Csums[i] > *Diameter){
			cout << "Diameter changed, before:" << *Diameter;
			*Diameter = Csums[i];
			cout << " after " << *Diameter << '\n';
		}
	}
}

Solstice_output* Solstice(OSwitchDemandMatrix* DemandM, uint64_t delay_ns, uint64_t circuitrate, uint64_t packetrate){
	uint64_t diameter = 0;//Fi = max of(Ri, Ci)
	Solstice_output* __res;
	/*calculate the duration nums(assume 0) and so on*/
	__res = new Solstice_output(used_ports, used_ports,0);

	//E <- D
	cout << "The packet demand is:\n";
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			__res->PSwitchMatrix[i][j] = DemandM->m_array[i][j];
			cout << __res->PSwitchMatrix[i][j] << ' ';
		}
		cout << '\n';
	}
    
	//D' <- D
	cout << "Trying Stuffing.\n";
	OSwitchDemandMatrix* StuffedDM = QuickStuff(DemandM);
    
	//T <- 0 in ns
	uint64_t total_time= 0;
	//r <- largest power of 2 smaller than max(D')
	OSwitchDemandMatrixElem* StuffMax = StuffedDM->Find_Max();

	uint64_t threshold = 1;
	while(threshold < StuffMax->value){
		threshold = threshold * 2;
	}
	threshold = threshold / 2;
	
	// i <- 0
	int dura_nums = 0;

	// while exsist row or column sum of D' > rpT do
	//caculate Ri and Ci
	SolsticeUpdateDiameter(StuffedDM, &diameter);
	
	//1Gbps = 1bpns
	uint64_t packetrate_bpns = packetrate;
	uint64_t circuitrate_bpns = circuitrate;
	
	while (diameter > packetrate_bpns * total_time){
		cout << "diameter is " << diameter << '\n';
		SliceOutput* slicetemp = BigSlice(StuffedDM, threshold);
		if(slicetemp->valid){
			cout << "Slicing is valid. threshold is" << threshold <<"\n";
			//ti <- min{D'(a,b) | Pi(a,b) = 1} / rc
			
			uint64_t min_dura_bytes = 0xffffffffffffffff;
			for (int i = 0; i < used_ports; i++){
				for (int j = 0; j < used_ports; j++){
					if(slicetemp->Permutation_Matrix[i][j]){
						min_dura_bytes = std::min(StuffedDM->m_array[i][j], min_dura_bytes);
					}
				}
			}
			cout << "min_dura_bytes is" << min_dura_bytes << '\n';
			uint64_t dura_time = min_dura_bytes / circuitrate_bpns;
			cout << "dura_time is " << dura_time << '\n';
			cout << "circuitrate_bpns is " << circuitrate_bpns << '\n';
			//__res->duration_time[dura_nums] = dura_time;
			__res->duration_time.push_back(dura_time);


			//D' <- D' - rc*ti*Pi  E <- E - rc*ti*Pi E <- ZeroEntriesBelow(E,0)
			for (int i = 0; i < used_ports; i++){
				for (int j = 0; j < used_ports; j++){
					if(slicetemp->Permutation_Matrix[i][j]){
						StuffedDM->m_array[i][j] = StuffedDM->m_array[i][j] - circuitrate_bpns * __res->duration_time[dura_nums];
						if (__res->PSwitchMatrix[i][j] < circuitrate_bpns * __res->duration_time[dura_nums]){
							__res->PSwitchMatrix[i][j] = 0;
						}
                        else{
                            __res->PSwitchMatrix[i][j] = __res->PSwitchMatrix[i][j] - circuitrate_bpns * __res->duration_time[dura_nums];
                        }
					}
				}
			}
            cout << "The Packet demand matrix should be:\n";
            for (int i = 0; i < used_ports; i++){
				for (int j = 0; j < used_ports; j++){
                    cout << __res->PSwitchMatrix[i][j] << ' ';
				}
                cout << '\n';
			}
			// T <- T + ti + delta
			total_time = total_time + __res->duration_time[dura_nums] + delay_ns;
			cout << "total_time is" << total_time << '\n';
			// i <- i + 1
			dura_nums = dura_nums + 1;
			cout << "dura_nums is" << dura_nums << '\n';

			cout << "The stuffed demand matrix now is:\n";
			for (int i = 0; i < used_ports; i++){
				for (int j = 0; j < used_ports; j++){
					cout << StuffedDM->m_array[i][j] << ' ';
				}
				cout << '\n';
			}
			//sleep(5);
		}
		else{
			//r <- r / 2
			if(threshold > 1){
				threshold = threshold / 2;
			}
			cout << "Slicing is not valid. Threshold is" << threshold << '\n';
			//sleep(5);
		}
		SolsticeUpdateDiameter(StuffedDM, &diameter);
		
	}
	
	// m <- i and prepare to set the circuit configures.
	__res->duration_nums = dura_nums;
	__res->circuit_configures = new bool** [__res->duration_nums];
	for (int i = 0; i < __res->duration_nums; i++){
		__res->circuit_configures[i] = new bool*[used_ports];
	}
	for (int i = 0; i < __res->duration_nums; i++){
		for (int j = 0; j < used_ports; j++){
			__res->circuit_configures[i][j] = new bool[used_ports];
		}
	}
	cout << "__res's duration nums is " << __res->duration_nums << '\n';
	//Now fill the circuit configurations and durations, calculate again.
    
	cout << "----------------Second calculating-----------------------------------------------\n";
	diameter = 0;
	total_time= 0;
	OSwitchDemandMatrix* cStuffedDM = QuickStuff(DemandM);
	// for (int i = 0; i < used_ports; i++){
	// 	for (int j = 0; j < used_ports; j++){
	// 		cout << cStuffedDM->m_array[i][j] << ' ';
	// 	}
	// 	cout << '\n';
	// }
	//cout << "Trying to find StuffMax.\n";
	StuffMax = cStuffedDM->Find_Max();
	//cout << "StuffMax's value is" << StuffMax->value << '\n';
	threshold = 1;
	while(threshold < StuffMax->value){
		threshold = threshold * 2;
	}
	threshold = threshold / 2;
	cout << "The second calculate's start threshold is " << threshold << '\n';
	dura_nums = 0;
	SolsticeUpdateDiameter(cStuffedDM, &diameter);
	while(diameter > packetrate_bpns * total_time){
		SliceOutput* slicetemp = BigSlice(cStuffedDM, threshold);
		if(slicetemp->valid){
			uint64_t min_dura_bytes = 0xffffffffffffffff;
			for (int i = 0; i < used_ports; i++){
				for (int j = 0; j < used_ports; j++){
					if(slicetemp->Permutation_Matrix[i][j]){
						min_dura_bytes = std::min(cStuffedDM->m_array[i][j], min_dura_bytes);
					}
				}
			}
			uint64_t dura_time = min_dura_bytes / circuitrate_bpns;
			// Fill the output's circuit configures.
			for (int i = 0; i < used_ports; i++){
				for (int j = 0; j < used_ports; j++){
					__res->circuit_configures[dura_nums][i][j] = slicetemp->Permutation_Matrix[i][j];
				}
			}
			for (int i = 0; i < used_ports; i++){
				for (int j = 0; j < used_ports; j++){
					if(slicetemp->Permutation_Matrix[i][j]){
						cStuffedDM->m_array[i][j] = cStuffedDM->m_array[i][j] - circuitrate_bpns * __res->duration_time[dura_nums];
						// __res->PSwitchMatrix[i][j] = __res->PSwitchMatrix[i][j] - circuitrate_bpns * __res->duration_time[dura_nums];
						// if (__res->PSwitchMatrix[i][j] < 0){
						// 	__res->PSwitchMatrix[i][j] = 0;
						// }
					}
				}
			}
			total_time = total_time + __res->duration_time[dura_nums] + delay_ns;
			dura_nums = dura_nums + 1;
		}
		else{
			if(threshold > 1){
				threshold = threshold / 2;
			}
		}
		SolsticeUpdateDiameter(cStuffedDM, &diameter);
	}
	return __res;
}

int main(int argc, char const *argv[]){
	uint64_t delay;
	uint64_t prate, orate;
	cout << "Input the configure delay:";
	cin >> delay;
	cout << "Input the rate of packet switch route:";
	cin >> prate;
	cout << "Input the rate of optical switch route:";
	cin >> orate;
    cout << "Input the matrix's width:";
    cin >> used_ports;
    uint64_t** mtx;
    mtx = new uint64_t* [used_ports];
    for (int i = 0; i < used_ports; i++){
        mtx[i] = new uint64_t[used_ports];
    }

    cout << "Input the whole matrix:\n";

    for(int i = 0;i < used_ports; i++){
        for (int j = 0; j < used_ports; j++){
            cin >> mtx[i][j];
	    //mtx[i][j] *= 10;
        }
    }

    OSwitchDemandMatrix* inmtx;
    inmtx = new OSwitchDemandMatrix(used_ports, used_ports);
    inmtx->Set_matrix(mtx, used_ports, used_ports);
    cout << "The demand matrix:\n";
    for(int i = 0;i < used_ports; i++){
        for (int j = 0; j < used_ports; j++){
			cout << mtx[i][j] << ' ';
        }
		cout << '\n';
    }

	cout << "Trying Solsitce.\n";
	Solstice_output* my_result = Solstice(inmtx, delay, orate, prate);
	cout << "Solsitce Over.\n";
	my_result->display();
    return 0;
}
