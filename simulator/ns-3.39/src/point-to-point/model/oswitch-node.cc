#include "ns3/ipv4.h"
#include "ns3/packet.h"
#include "ns3/ipv4-header.h"
#include "ns3/pause-header.h"
#include "ns3/interface-tag.h"
#include "ns3/boolean.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "oswitch-node.h"
#include "qbb-net-device.h"
#include "ppp-header.h"
#include "ns3/int-header.h"
#include "ns3/simulator.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstring>
#include <iostream>
#include "ns3/tcp-header.h"
#include "ns3/udp-header.h"
#include "ns3/custom-priority-tag.h"
#include "ns3/feedback-tag.h"
#include "ns3/unsched-tag.h"

namespace ns3 {

TypeId OSwitchNode::GetTypeId (void)
{
    static TypeId tid = TypeId ("ns3::OSwitchNode")
	                    .SetParent<Node> ()
	                    .AddConstructor<OSwitchNode> ()
						.AddAttribute("ReconfigureDelay",
									  "Reconfig delay time in ns.",
									  UintegerValue(5000000),
									  MakeUintegerAccessor(&OSwitchNode::m_reconf_delay),
									  MakeUintegerChecker<uint64_t>())
						.AddAttribute("OswitchRate",
									  "Optical switch's rate.",
									  UintegerValue(400000000000),
									  MakeUintegerAccessor(&OSwitchNode::m_rate),
									  MakeUintegerChecker<uint64_t>())
						.AddAttribute("PswitchRate",
									  "Packet switch's rate.",
									  UintegerValue(25000000000),
									  MakeUintegerAccessor(&OSwitchNode::p_rate),
									  MakeUintegerChecker<uint64_t>())
						.AddAttribute("DurationExtentionRate",
									  "Duration time extention Rate (theta).",
									  DoubleValue(1.5),
									  MakeDoubleAccessor(&OSwitchNode::theta),
									  MakeDoubleChecker<double>())
	                    ;
    return tid;
}

OSwitchNode::OSwitchNode() {
    m_node_type = 3;
    for (uint32_t i = 0; i < pCnt; i++)
		for (uint32_t j = 0; j < pCnt; j++)
			m_bytes[i][j] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_txBytes[i] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_lastPktSize[i] = m_lastPktTs[i] = 0;
	for (uint32_t i = 0; i < pCnt; i++)
		m_u[i] = 0;
	m_DemandM = new OSwitchDemandMatrix;
	for (uint32_t i = 0; i < pCnt; i++){
		for (uint32_t j = 0; j < pCnt; j++){
			m_circuits[i][j] = 0;
			m_circuits[j][i] = 0;
		}
	}
}

int OSwitchNode::GetOutDev(Ptr<const Packet> p, CustomHeader &ch) {
	// look up entries
	Ptr<Packet> cp = p->Copy();

	PppHeader ph; cp->RemoveHeader(ph);
	Ipv4Header ih; cp->RemoveHeader(ih);
	auto entry = m_rtTable.find(ih.GetDestination().Get());

	// no matching entry
	if (entry == m_rtTable.end())
		return -1;
	// entry found
	auto &nexthops = entry->second;

	return nexthops[0];

}

int OSwitchNode::GetInDev(Ptr<const Packet> p, CustomHeader &ch) {
	Ptr<Packet> cp = p->Copy();

	PppHeader ph; cp->RemoveHeader(ph);
	Ipv4Header ih; cp->RemoveHeader(ih);
	auto entry = m_rtTable.find(ih.GetSource().Get());

	// no matching entry
	if (entry == m_rtTable.end())
		return -1;
	// entry found
	auto &lasthop = entry->second;

	return lasthop[0];
}


void OSwitchNode::SendToDev(Ptr<Packet> p, CustomHeader &ch){
	int in_idx = GetInDev(p, ch);
	int out_idx = GetOutDev(p, ch);

	//std::cout << "At time " << Simulator::Now().GetNanoSeconds() << "Optical Switch trying to send packet.\n";

	if (out_idx >= 0 && in_idx >= 0) {
		NS_ASSERT_MSG(m_devices[out_idx]->IsLinkUp(), "The routing table look up should return link that is up");

		if (m_circuits[in_idx-1][out_idx-1] == 1 /*&& m_circuits[out_idx-1][in_idx-1] == 1*/)
		{	// The circuit is connected

			// determine the qIndex
			uint32_t qIndex=1;
			m_bytes[in_idx][out_idx] += p->GetSize();
			//std::cout << "Optical Switch " << this->GetId() << " is sending via port:" << in_idx <<" to port:" << out_idx << '\n';
			m_devices[out_idx]->SwitchSend(qIndex, p, ch);
			DynamicCast<QbbNetDevice>(m_devices[out_idx])->totalBytesRcvd += p->GetSize(); // Attention: this is the egress port's total received packets. Not the ingress port.
		}
		//else std::cout << "Optical Switch " << this->GetId() << " 's src port:" << in_idx <<" to port:" << out_idx << "is not connected!" << '\n' ;
		{	//The circuit is not connected

			return; // Drop
		}
	} else
		std::cout << "outdev and indev not found! Dropped. This should not happen. Debugging required!" << std::endl;
	return; // Drop
}

void OSwitchNode::AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx){
	uint32_t dip = dstAddr.Get();
	m_rtTable[dip].push_back(intf_idx);
}

void OSwitchNode::ShowTableEntry(){
	std::cout << "The optical switch's Table Entry:\n";
	for (auto i = m_rtTable.begin(); i != m_rtTable.end(); i++){
		uint32_t dip = i->first;
		std::cout << "\tTo the des: " << dip << ":\n\t\tInterface:";
		for (auto j = i->second.begin(); j != i->second.end(); j++){
			std::cout << *j << ' ';
		}
		std::cout << '\n';
	}
}

void OSwitchNode::ClearTable() {
	m_rtTable.clear();
}

void OSwitchNode::CalculateConfigures(){
	m_conf_and_times = Solstice(m_DemandM, m_reconf_delay / 10000 , m_rate / 1000000000, p_rate / 1000000000);
	std::cout << "Solstice over.\n";
	DayorNight = 0;
	curr_dura_num = 0;
}

void OSwitchNode::Reconfigure() {
	if (DayorNight == 0){	// Night time. Do the reconfigure work.
		Time dur_time = SetNewTopology();
		DayorNight = 1;	//Turn to Daytime
		std::cout << "At Time:" << Simulator::Now().GetNanoSeconds();
		std::cout << " Night time to day time. Calling Simulator::Schedule.\n";
		Simulator::Schedule(dur_time, &OSwitchNode::Reconfigure, this);
	}
	else {	//Day time.	Shutdown evrey port and prepare to reconfigure.
		//Cut all the circuits
		ClearTopology();
		DayorNight = 0;
		Time reconfiguretime = Time::FromInteger(m_reconf_delay, Time::NS);
		std::cout << "At Time:" << Simulator::Now().GetNanoSeconds();
		std::cout << " Day time to night time. Calling Simulator::Schedule.\n";
		Simulator::Schedule(reconfiguretime, &OSwitchNode::Reconfigure, this);
	}
}


void OSwitchNode::SwitchNotifyDequeue(uint32_t ifIndex, uint32_t qIndex, Ptr<Packet> p) {
	m_txBytes[ifIndex] += p->GetSize();
	m_lastPktSize[ifIndex] = p->GetSize();
	m_lastPktTs[ifIndex] = Simulator::Now().GetTimeStep();
}

bool  OSwitchNode::SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch) {
	SendToDev(packet, ch);
	return true;
}

void OSwitchNode::ClearTopology(){
	for (int i = 0; i < pCnt; i++){
		for (int j = 0; j < pCnt; j++){
			m_circuits[i][j] = 0;
			m_circuits[j][i] = 0;
		}
	}
}

Time OSwitchNode::SetNewTopology(){
	std::cout<< "SetNewTopology Called.\n";
	Time _res;
	if(curr_dura_num >= m_conf_and_times->duration_nums){
		std::cout << "curr_dura_num exceeded duration_nums, just return.\n";
		uint64_t this_time = m_conf_and_times->duration_time[curr_dura_num - 1];
		_res = Time::FromInteger(this_time, Time::NS);
		return _res;
	}
	uint64_t this_time = m_conf_and_times->duration_time[curr_dura_num];
	if (theta != 1.0){
		this_time = (uint64_t)(this_time * theta);
	}
	//std::cout << "m_circuits is:\n";
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			if(m_conf_and_times->circuit_configures[curr_dura_num][i][j] == 1){
				m_circuits[i][j] = 1;
				m_circuits[j][i] = 1;
			}
			else{
				m_circuits[i][j] = 0;
				m_circuits[j][i] = 0;
			}
			//std::cout << m_circuits[i][j] << ' ';
		}
		//std::cout << '\n';
	}
	std::cout<< "duration time is: " << this_time <<".\n";
	_res = Time::FromInteger(this_time, Time::NS);
	curr_dura_num = curr_dura_num + 1;
	return _res;
}

void OSwitchNode::InitDemandM(int rnum, int cnum){
	std::cout << "Initing Demand Matrix, rnum:" << rnum << " cnum:" << cnum << '\n';
	m_DemandM->Numrow = rnum;
	m_DemandM->Numcol = cnum;
	m_DemandM->m_array = new uint64_t* [m_DemandM->Numrow];
	for (int i = 0;i < m_DemandM->Numrow; i++){
		m_DemandM->m_array[i] = new uint64_t[m_DemandM->Numcol];
	}
	for (int i = 0; i < m_DemandM->Numrow; i++){
		for (int j = 0; j < m_DemandM->Numcol; j++){
			m_DemandM->m_array[i][j] = 0;
		}
	}
}

void OSwitchNode::AddDemandM(int inport, int outport, uint64_t value){
	m_DemandM->m_array[inport-1][outport-1] += value;
}

void OSwitchNode::DemandMSetup(int used_ports_num){
	this->used_ports = used_ports_num;
	InitDemandM(used_ports_num, used_ports_num);
}

void OSwitchNode::SetFlowDemand(int in_indx, int out_indx, uint64_t size){
	AddDemandM(in_indx, out_indx, size);
}

void OSwitchNode::DisplayDemandM(){
	std::cout << "Optical Switch's Demand matrix:\n";
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			std::cout << m_DemandM->m_array[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

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


OSwitchDemandMatrix* OSwitchNode::QuickStuff(OSwitchDemandMatrix* matrix){
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


SliceOutput* OSwitchNode::BigSlice(OSwitchDemandMatrix* matrix, uint64_t threshold){
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

void OSwitchNode::SolsticeUpdateDiameter(OSwitchDemandMatrix* DemandM, uint64_t* Diameter){
	*Diameter = 0;
	uint64_t Rsums[used_ports];//Ri
	uint64_t Csums[used_ports];//Ci
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
			*Diameter = Rsums[i];
		}
		if(Csums[i] > *Diameter){
			*Diameter = Csums[i];
		}
	}
}

Solstice_output* OSwitchNode::Solstice(OSwitchDemandMatrix* DemandM, uint64_t delay_ns, uint64_t circuitrate, uint64_t packetrate){
	std::cout << "Solstice Called. " << "delay in ns is " << delay_ns << " circuitrate is " << circuitrate << " packetrate is " << packetrate << "\n";
	uint64_t diameter = 0;//Fi = max of(Ri, Ci)
	Solstice_output* __res;
	/*calculate the duration nums(assume 0) and so on*/
	__res = new Solstice_output(used_ports, used_ports,0);
	
	//E <- D
	for (int i = 0; i < used_ports; i++){
		for (int j = 0; j < used_ports; j++){
			__res->PSwitchMatrix[i][j] = DemandM->m_array[i][j];
		}
	}

	//D' <- D
	OSwitchDemandMatrix* StuffedDM = QuickStuff(DemandM);
	std::cout << "QuickStuff Over.\n";
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
		SliceOutput* slicetemp = BigSlice(StuffedDM, threshold);
		if(slicetemp->valid){
			//ti <- min{D'(a,b) | Pi(a,b) = 1} / rc
			uint64_t min_dura_bytes = 0xffffffffffffffff;
			for (int i = 0; i < used_ports; i++){
				for (int j = 0; j < used_ports; j++){
					if(slicetemp->Permutation_Matrix[i][j]){
						min_dura_bytes = std::min(StuffedDM->m_array[i][j], min_dura_bytes);
					}
				}
			}
			uint64_t dura_time = min_dura_bytes / circuitrate_bpns;
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
						else {
							__res->PSwitchMatrix[i][j] = __res->PSwitchMatrix[i][j] - circuitrate_bpns * __res->duration_time[dura_nums];
						}
					}
				}
			}
			// T <- T + ti + delta
			total_time = total_time + __res->duration_time[dura_nums] + delay_ns;
			// i <- i + 1
			dura_nums = dura_nums + 1;
		}
		else{
			//r <- r / 2
			if (threshold > 1){
				threshold = threshold / 2;
			}
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
	std::cout << "Solstice First Calculation Over.\n";
	diameter = 0;
	total_time= 0;
	OSwitchDemandMatrix* cStuffedDM = QuickStuff(DemandM);
	StuffMax = cStuffedDM->Find_Max();
	threshold = 1;
	while(threshold < StuffMax->value){
		threshold = threshold * 2;
	}
	threshold = threshold / 2;
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
bool OSwitchNode::CircuitConnected(Ptr<Packet>p, CustomHeader &ch) {
	int in_idx = GetInDev(p, ch);
	int out_idx = GetOutDev(p, ch);
	//std::cout << "At Time:" << Simulator::Now().GetNanoSeconds();
	//std::cout << " Optical's CircuitConnected is called. in_idx:" << in_idx << " out_idx:" << out_idx << '\n';
	bool result = false;
	//std::cout << "m_circuits["<< in_idx-1 << "]["<<out_idx-1<<"]" << "is" <<m_circuits[in_idx-1][out_idx-1];
	//std::cout << "m_circuits["<<out_idx-1<<"]["<<in_idx-1<<"]" << "is" <<m_circuits[out_idx-1][in_idx-1];
	if (m_circuits[in_idx-1][out_idx-1] == 1 && m_circuits[out_idx-1][in_idx-1] == 1){
		result = true;
	}
	//std::cout << "result is " << result << '\n';
	return result;
}

void OSwitchNode::man_init_configures(int tornum, int duration_nums){
	m_conf_and_times = new Solstice_output(tornum, tornum, duration_nums);
	m_conf_and_times->circuit_configures = new bool** [m_conf_and_times->duration_nums];
	for (int i = 0; i < m_conf_and_times->duration_nums; i++){
		m_conf_and_times->circuit_configures[i] = new bool*[tornum];
	}
	for (int i = 0; i < m_conf_and_times->duration_nums; i++){
		for (int j = 0; j < tornum; j++){
			m_conf_and_times->circuit_configures[i][j] = new bool[tornum];
		}
	}
}

void OSwitchNode::man_insert_duration_time(uint64_t dura_time) {
	m_conf_and_times->duration_time.push_back(dura_time);
}

void OSwitchNode::man_set_configures( int curr_num, int row, int col, bool value){
	m_conf_and_times->circuit_configures[curr_num][row][col] = value;
	m_conf_and_times->circuit_configures[curr_num][col][row] = value;
}

}/*namespace ns3*/