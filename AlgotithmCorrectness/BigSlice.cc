#include <iostream>
#include <vector>
#include <cstring>
#include <algorithm>

using namespace std;

int used_ports;


struct OSwitchDemandMatrixElem {
	int row;
	int col;
	int value;
};

class OSwitchDemandMatrix {
public:
	int Numrow;
	int Numcol;
	int** m_array;

	OSwitchDemandMatrix(int row, int col){
		Numrow = row;
		Numcol = col;
		m_array = new int* [Numrow];
		for (int i = 0;i < Numrow; i++){
			m_array[i] = new int[Numcol];
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

	OSwitchDemandMatrixElem Find_Max(){
		int _max_value = 0;
		OSwitchDemandMatrixElem __res;
		__res.col = __res.row = 0;
		__res.value = 0;
		for (int i = 0; i < Numrow; i++)
		{
			for (int j = 0; j < Numcol; j++)
			{
				if (m_array[i][j] > _max_value){
					__res.row = i;
					__res.col = j;
					__res.value = m_array[i][j];
					_max_value = m_array[i][j];
				}
			}
		}
        return __res;
	}

	void Set_matrix(int** mtx, int rnum, int cnum){
		for(int i = 0;i < rnum; i++){
			for (int j = 0;j < cnum; j++){
				m_array[i][j] = mtx[i][j];
			}
		}
	}

    void Print_matrix(){
        for(int i = 0;i < Numrow; i++){
			for (int j = 0;j < Numcol; j++){
				cout << m_array[i][j] << ' ';
			}
            cout << '\n';
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
    void print() {
        if(valid){
            cout << "This is a valid result.\n";
        }
        for (int i = 0; i < rowNum; i++){
            for (int j = 0; j < colNum; j++){
                cout << Permutation_Matrix[i][j] << ' ';
            }
            cout << '\n';
        }
    }
};

SliceOutput* BigSlice(OSwitchDemandMatrix* matrix, int threshold){
	SliceOutput* __res;
	//Return value, default to be not valid, determine it will be valid or not later.
	__res = new SliceOutput(false, used_ports, used_ports);

	//Copy the matrix
	int** tempmatrix = new int*[used_ports];
	for (int i = 0;i < used_ports; i++){
		tempmatrix[i] = new int[used_ports];
	}
	for (int i = 0;i < used_ports; i++){
		for (int j = 0;j < used_ports; j++){
			tempmatrix[i][j] = matrix->m_array[i][j];
		}
	}
//-------------------------------------------------------------------------
	//First Do D'' <- ZeroEntiresBelow(D', r) AND B-<BinaryMatrixOf(D'')
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			if (tempmatrix[i][j] < threshold){
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
    cout << "result sum is: " << resultsum << '\n';
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
	int Rsums[used_ports];//Ri
	int Csums[used_ports];//Ci
	int diameter;//Fi = max of(Ri, Ci)
	OSwitchDemandMatrix* __res;
	__res = new OSwitchDemandMatrix(used_ports, used_ports);
	int** tempmatrix = new int*[used_ports];
	for (int i = 0;i < used_ports; i++){
		tempmatrix[i] = new int[used_ports];
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
	
	int temp;
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

int main(int argc, char const *argv[]){
    cout << "Input the matrix's width:";
    cin >> used_ports;
    int** mtx;
    mtx = new int* [used_ports];
    for (int i = 0; i < used_ports; i++){
        mtx[i] = new int[used_ports];
    }

    cout << "Input the whole matrix:\n";

    for(int i = 0;i < used_ports; i++){
        for (int j = 0; j < used_ports; j++){
            cin >> mtx[i][j];
        }
    }
    OSwitchDemandMatrix* inmtx;
    inmtx = new OSwitchDemandMatrix(used_ports, used_ports);
    inmtx->Set_matrix(mtx, used_ports, used_ports);
    cout << "The matrix before stuff:\n";
    inmtx->Print_matrix();

    OSwitchDemandMatrix* outmtx = QuickStuff(inmtx);
    cout << "The matrix after stuff:\n";
    outmtx->Print_matrix();

    cout << "Now try to do the Big Slice. r=64.\n";
    SliceOutput* result;
    result = BigSlice(outmtx, 64);
    result->print();
    return 0;
}
