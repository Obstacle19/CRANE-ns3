#ifndef OSWITCH_NODE_H
#define OSWITCH_NODE_H

#include <iostream>
#include <unordered_map>
#include <vector>
#include <ns3/node.h>
#include "qbb-net-device.h"
#include "switch-mmu.h"
#include "pint.h"

namespace ns3 {

class Packet;

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
	OSwitchDemandMatrix(){}
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

	void Set_circuit_configures(uint64_t*** circuit_conf, int dura_times, int rnum, int cnum){
		for (int dt = 0;dt < dura_times; dt++){
			for(int i = 0;i < rnum; i++){
				for (int j = 0;j < cnum; j++){
					circuit_configures[dt][i][j] = circuit_conf[dt][i][j];
				}
			}
		}
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
		std::cout << "There are " << duration_nums << " duration nums:\n";
		for (int i = 0; i < duration_nums; i++){
			std::cout << "The number " << i << " st/nd/th duation is:\n";
			for (int a = 0; a < Numrow; a++){
				for(int b = 0; b < Numcol; b++){
					std::cout<< circuit_configures[i][a][b] << ' ';
				}
				std::cout << '\n';
			}
			std::cout << "This configure's duration is " << duration_time[i] << " ns.\n";
		}
		
		std::cout << "And the Packet demand matrix is:\n";
		for (int a = 0; a < Numrow; a++){
			for(int b = 0; b < Numcol; b++){
				std::cout<< PSwitchMatrix[a][b] << ' ';
			}
			std::cout << '\n';
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


class OSwitchNode : public Node{
    static const uint32_t pCnt = 257;	// Number of ports used
    std::unordered_map<uint32_t, std::vector<int> > m_rtTable; // map from ip address (u32) to possible ECMP port (index of dev)

    uint32_t m_bytes[pCnt][pCnt]; // m_bytes[inDev][outDev][qidx] is the bytes from inDev enqueued for outDev at qidx
	
	uint64_t m_txBytes[pCnt]; // counter of tx bytes
	uint32_t m_circuits[pCnt][pCnt];

	uint32_t m_lastPktSize[pCnt];
	uint64_t m_lastPktTs[pCnt]; // ns
	double m_u[pCnt];
	double theta;

	OSwitchDemandMatrix* m_DemandM;
	int used_ports;		//How many ports uesd.
	int curr_dura_num;	//The dura_num of now. from 0 to m_conf_and_times's max.
	//Time m_daytime;		//Each Day's time.
	//Time m_reconfiguretime;	
	Solstice_output* m_conf_and_times;	//Solstice algorithm's output.
	uint64_t m_reconf_delay;			//reconfigure's delay.
	uint64_t m_rate;					//optical switch's rate.
	uint64_t p_rate;					//packet switch's rate.
	bool DayorNight; // 1 for day and 0 for night.

private:
	int GetOutDev(Ptr<const Packet>, CustomHeader &ch);
	int GetInDev(Ptr<const Packet>, CustomHeader &ch);
	void SendToDev(Ptr<Packet>p, CustomHeader &ch);
	
public:
	static TypeId GetTypeId (void);
	OSwitchNode();
	void AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx);
	void ShowTableEntry();
	void ClearTable();
	void DemandMSetup(int used_ports_num);
	void SetFlowDemand(int in_indx, int out_indx, uint64_t size);
	void DisplayDemandM();
    void Reconfigure();
	bool SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch);
	void SwitchNotifyDequeue(uint32_t ifIndex, uint32_t qIndex, Ptr<Packet> p);
	void ClearTopology();
	OSwitchDemandMatrix* QuickStuff(OSwitchDemandMatrix* matrix);
	SliceOutput* BigSlice(OSwitchDemandMatrix* matrix, uint64_t threshold);
	void SolsticeUpdateDiameter(OSwitchDemandMatrix* DemandM, uint64_t* Diameter);
	void CalculateConfigures();
	Solstice_output* Solstice(OSwitchDemandMatrix* DemandM, uint64_t delay_ns, uint64_t circuitrate_Gps, uint64_t packetrate_Gps);
	Time SetNewTopology();
	void InitDemandM(int rnum, int cnum);
	void AddDemandM(int inport, int outport, uint64_t value);
	bool CircuitConnected(Ptr<Packet>p, CustomHeader &ch);
	void display_configures(){
		m_conf_and_times->display();
	}
	void SetDayorNight(bool value) {
		DayorNight = value;
	}
	void man_init_configures(int tornum, int duration_nums);
	void man_insert_duration_time(uint64_t dura_time);
	void man_set_configures(int curr_num, int row, int col, bool value);
};

}/*namespace ns3*/

#endif /*OSWITCH_NODE_H*/