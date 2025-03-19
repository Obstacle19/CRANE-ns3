#include <iostream>
#include <algorithm>
using namespace std;

int used_ports;

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
				temp = diameter - max(Rsums[i], Csums[j]);
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
    cout << "Input the matrix's width.:";
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
    return 0;
}
