/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License version 2 as
* published by the Free Software Foundation;
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#undef PGO_TRAINING
#define PATH_TO_PGO_CONFIG "path_to_pgo_config"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <time.h>
#include "ns3/core-module.h"
#include "ns3/qbb-helper.h"
#include "ns3/point-to-point-helper.h"
#include "ns3/applications-module.h"
#include "ns3/internet-module.h"
#include "ns3/global-route-manager.h"
#include "ns3/ipv4-static-routing-helper.h"
#include "ns3/packet.h"
#include "ns3/error-model.h"
#include <ns3/rdma.h>
#include <ns3/rdma-client.h>
#include <ns3/rdma-client-helper.h>
#include <ns3/rdma-driver.h>
#include <ns3/switch-node.h>
#include <ns3/sim-setting.h>
#include <ns3/oswitch-node.h>
#include <ns3/torswitch-node.h>

using namespace ns3;
using namespace std;

NS_LOG_COMPONENT_DEFINE("GENERIC_SIMULATION");

uint32_t cc_mode = 1; // congestion control ，用于指定拥塞控制模式的选择，默认为 1
bool enable_qcn = true; // 指示是否启用 QCN
// 数据包的有效负载大小、链路层数据块的大小、链路层确认（ACK）发送的时间间隔
uint32_t packet_payload_size = 1000, l2_chunk_size = 0, l2_ack_interval = 0;
// 模拟中暂停的时间、模拟运行的结束时间
double pause_time = 5, simulator_stop_time = 3.01;
// 分别存储数据速率、链路延迟、拓扑文件、流文件、跟踪文件和跟踪输出文件的路径或名称
std::string data_rate, link_delay, topology_file, flow_file, trace_file, trace_output_file;
// 分别用于存储完成时间（FCT）和优先流控制（PFC）输出结果的文件名
std::string fct_output_file = "fct.txt";
std::string pfc_output_file = "pfc.txt";

double alpha_resume_interval = 55, rp_timer, ewma_gain = 1 / 16;
double rate_decrease_interval = 4;
uint32_t fast_recovery_times = 5;
std::string rate_ai, rate_hai, min_rate = "100Mb/s";
std::string dctcp_rate_ai = "1000Mb/s";

bool clamp_target_rate = false, l2_back_to_zero = false;
double error_rate_per_link = 0.0;
uint32_t has_win = 1;
uint32_t global_t = 1;
uint32_t mi_thresh = 5;
bool var_win = false, fast_react = true;
bool multi_rate = true;
bool sample_feedback = false;
double pint_log_base = 1.05;
double pint_prob = 1.0;
double u_target = 0.95;
uint32_t int_multi = 1;
bool rate_bound = true;

uint32_t ack_high_prio = 0;
uint64_t link_down_time = 0;
uint32_t link_down_A = 0, link_down_B = 0;

uint32_t enable_trace = 1;

uint32_t buffer_size = 16;

uint32_t qlen_dump_interval = 100000000, qlen_mon_interval = 100;
uint64_t qlen_mon_start = 2000000000, qlen_mon_end = 2100000000;
string qlen_mon_file;
string optical_conf_file;

unordered_map<uint64_t, uint32_t> rate2kmax, rate2kmin;
unordered_map<uint64_t, double> rate2pmax;

uint64_t optical_reconfigure_delay; // OCS 的重配置延迟
uint64_t optical_oswitch_rate; // 还是 OCS 的传输速率，只是复制了一遍，用于 Solstice 算法的输入
uint64_t optical_pswitch_rate; // 还是 EPS 的传输速率，只是复制了一遍，用于 Solstice 算法的输入
double optical_theta;

#define MAN_OPTICAL 1

/************************************************
 * Runtime varibles
 ***********************************************/

std::ifstream topof, flowf, tracef; // 文件流变量，用于读取不同的配置文件

NodeContainer n; // 一个容器，用于存储模拟中的所有节点
NetDeviceContainer switchToSwitchInterfaces;
std::map< uint32_t, std::map< uint32_t, std::vector<Ptr<QbbNetDevice>> > > switchToSwitch;

std::map<uint32_t, uint32_t> switchNumToId;
std::map<uint32_t, uint32_t> switchIdToNum;
std::map<uint32_t, NetDeviceContainer> switchUp;
std::map<uint32_t, NetDeviceContainer> switchDown;
//NetDeviceContainer switchUp[switch_num];
std::map<uint32_t, NetDeviceContainer> sourceNodes;
std::map<Ptr<Node>, std::vector<Ptr<Node> > > Tor2Server;

// NodeContainer servers;
// NodeContainer tors;

uint64_t nic_rate;

uint64_t maxRtt, maxBdp;

struct Interface {
	uint32_t idx;
	bool up;
	uint64_t delay;
	uint64_t bw;

	Interface() : idx(0), up(false) {}
};

map<Ptr<Node>, map<Ptr<Node>, Interface> > nbr2if;
// Mapping destination to next hop for each node: <node, <dest, <nexthop0, ...> > >
map<Ptr<Node>, map<Ptr<Node>, vector<Ptr<Node> > > > nextHop;
map<Ptr<Node>, map<Ptr<Node>, uint64_t> > pairDelay;
map<Ptr<Node>, map<Ptr<Node>, uint64_t> > pairTxDelay;
map<uint32_t, map<uint32_t, uint64_t> > pairBw;
map<Ptr<Node>, map<Ptr<Node>, uint64_t> > pairBdp;
map<uint32_t, map<uint32_t, uint64_t> > pairRtt;

std::vector<Ipv4Address> serverAddress;

std::unordered_map<uint32_t, unordered_map<uint32_t, uint16_t> > portNumder;

struct FlowInput {
	uint64_t src, dst, pg, maxPacketCount, port, dport;
	double start_time;
	uint32_t idx;
};

FlowInput flow_input = {0};
uint32_t flow_num;

Ptr<OSwitchNode> OpticalSwitch;

void ReadFlowInput() {
	cout << "Going to read flow input.\n";
	if (flow_input.idx < flow_num) {
		flowf >> flow_input.src >> flow_input.dst >> flow_input.pg >> flow_input.dport >> flow_input.maxPacketCount >> flow_input.start_time;
		std::cout << "Flow " << flow_input.src << " " << flow_input.dst << " " << flow_input.pg << " " << flow_input.dport << " " << flow_input.maxPacketCount << " " << flow_input.start_time << " " << Simulator::Now().GetSeconds() << std::endl;
		NS_ASSERT(n.Get(flow_input.src)->GetNodeType() == 0 && n.Get(flow_input.dst)->GetNodeType() == 0);
		//cout << "ASSERT over.\n";
	}
}

void ScheduleFlowInputs() {
	cout << "Called ScheduleFlowInputs().\n";
	while (flow_input.idx < flow_num && Seconds(flow_input.start_time) <= Simulator::Now()) {
		uint32_t port = portNumder[flow_input.src][flow_input.dst]++; // get a new port number
		RdmaClientHelper clientHelper(flow_input.pg, serverAddress[flow_input.src], serverAddress[flow_input.dst], port, flow_input.dport, flow_input.maxPacketCount, has_win ? (global_t == 1 ? maxBdp : pairBdp[n.Get(flow_input.src)][n.Get(flow_input.dst)]) : 0, global_t == 1 ? maxRtt : pairRtt[flow_input.src][flow_input.dst], Simulator::GetMaximumSimulationTime());
		//cout << "appCon is installing in " << flow_input.src << '\n';
		ApplicationContainer appCon = clientHelper.Install(n.Get(flow_input.src));
		//cout << "appCon is installed in " << n.Get(flow_input.src) << '\n';
//		appCon.Start(Seconds(flow_input.start_time));
		appCon.Start(Seconds(0)); // setting the correct time here conflicts with Sim time since there is already a schedule event that triggered this function at desired time.
		// get the next flow input
		flow_input.idx++;
		ReadFlowInput();
	}

	// schedule the next time to run this function
	if (flow_input.idx < flow_num) {
		cout << "schedule the next time to run ScheduleFlowInputs.\n";
		Simulator::Schedule(Seconds(flow_input.start_time) - Simulator::Now(), ScheduleFlowInputs);
	} else { // no more flows, close the file
		flowf.close();
	}
}

Ipv4Address node_id_to_ip(uint32_t id) {
	return Ipv4Address(0x0b000001 + ((id / 256) * 0x00010000) + ((id % 256) * 0x00000100));
}

uint32_t ip_to_node_id(Ipv4Address ip) {
	return (ip.Get() >> 8) & 0xffff;
}

void qp_finish(FILE* fout, Ptr<RdmaQueuePair> q) {
	uint32_t sid = ip_to_node_id(q->sip), did = ip_to_node_id(q->dip);
	uint64_t base_rtt = pairRtt[sid][did], b = pairBw[sid][did];
	uint32_t total_bytes = q->m_size + ((q->m_size - 1) / packet_payload_size + 1) * (CustomHeader::GetStaticWholeHeaderSize() - IntHeader::GetStaticSize()); // translate to the minimum bytes required (with header but no INT)
	uint64_t standalone_fct = base_rtt + total_bytes * 8000000000lu / b;
	// sip, dip, sport, dport, size (B), start_time, fct (ns), standalone_fct (ns)
	fprintf(fout, "%08x %08x %u %u %lu %lu %lu %lu\n", q->sip.Get(), q->dip.Get(), q->sport, q->dport, q->m_size, q->startTime.GetTimeStep(), (Simulator::Now() - q->startTime).GetTimeStep(), standalone_fct);
	fflush(fout);

	// remove rxQp from the receiver
	Ptr<Node> dstNode = n.Get(did);
	Ptr<RdmaDriver> rdma = dstNode->GetObject<RdmaDriver> ();
	rdma->m_rdma->DeleteRxQp(q->sip.Get(), q->m_pg, q->sport);
}

void get_pfc(FILE* fout, Ptr<QbbNetDevice> dev, uint32_t type) {
	fprintf(fout, "%lu %u %u %u %u\n", Simulator::Now().GetTimeStep(), dev->GetNode()->GetId(), dev->GetNode()->GetNodeType(), dev->GetIfIndex(), type);
}

struct QlenDistribution {
	vector<uint32_t> cnt; // cnt[i] is the number of times that the queue len is i KB

	void add(uint32_t qlen) {
		uint32_t kb = qlen / 1000;
		if (cnt.size() < kb + 1)
			cnt.resize(kb + 1);
		cnt[kb]++;
	}
};

map<uint32_t, map<uint32_t, QlenDistribution> > queue_result;

void monitor_buffer(FILE* qlen_output, NodeContainer *n, Ptr<OSwitchNode> mocs) {	//监视所有switch的buffer情况
	for (uint32_t i = 0; i < n->GetN(); i++) {
		if (n->Get(i)->GetNodeType() == 1 && n->Get(i) != mocs) { // is switch and not optical switch
			Ptr<SwitchNode> sw = DynamicCast<SwitchNode>(n->Get(i));
			if (queue_result.find(i) == queue_result.end())
				queue_result[i];
			for (uint32_t j = 1; j < sw->GetNDevices(); j++) {
				uint32_t size = 0;
				for (uint32_t k = 0; k < SwitchMmu::qCnt; k++)
					size += sw->m_mmu->egress_bytes[j][k];
				queue_result[i][j].add(size);
			}
		}
	}
	if (Simulator::Now().GetTimeStep() % qlen_dump_interval == 0) {
		fprintf(qlen_output, "time: %lu\n", Simulator::Now().GetTimeStep());
		for (auto &it0 : queue_result)
			for (auto &it1 : it0.second) {
				fprintf(qlen_output, "%u %u", it0.first, it1.first);
				auto &dist = it1.second.cnt;
				for (uint32_t i = 0; i < dist.size(); i++)
					fprintf(qlen_output, " %u", dist[i]);
				fprintf(qlen_output, "\n");
			}
		fflush(qlen_output);
	}
	if (Simulator::Now().GetTimeStep() < qlen_mon_end)
		Simulator::Schedule(NanoSeconds(qlen_mon_interval), &monitor_buffer, qlen_output, n, mocs);
}

void CalculateRoute(Ptr<Node> host) {	//使用广度优先寻找最短路径
	// queue for the BFS.
	vector<Ptr<Node> > q;
	// Distance from the host to each node.
	map<Ptr<Node>, int> dis;
	map<Ptr<Node>, uint64_t> delay;
	map<Ptr<Node>, uint64_t> txDelay;
	map<Ptr<Node>, uint64_t> bw;
	// init BFS.
	q.push_back(host);
	dis[host] = 0;
	delay[host] = 0;
	txDelay[host] = 0;
	bw[host] = 0xfffffffffffffffflu;
	// BFS.
	for (int i = 0; i < (int)q.size(); i++) {
		Ptr<Node> now = q[i];
		int d = dis[now];
		for (auto it = nbr2if[now].begin(); it != nbr2if[now].end(); it++) {
			// skip down link
			if (!it->second.up)
				continue;
			if (now->GetId() == 0){ // Avoid Optical switch
				continue;
			}
			Ptr<Node> next = it->first;
			// If 'next' have not been visited.
			if (dis.find(next) == dis.end()) {
				dis[next] = d + 1;
				delay[next] = delay[now] + it->second.delay;
				txDelay[next] = txDelay[now] + packet_payload_size * 1000000000lu * 8 / it->second.bw;
				bw[next] = std::min(bw[now], it->second.bw);
				// we only enqueue switch, because we do not want packets to go through host as middle point
				if (next->GetNodeType())
					q.push_back(next);
			}
			// if 'now' is on the shortest path from 'next' to 'host'.
			if (d + 1 == dis[next]) {
				nextHop[next][host].push_back(now);
			}
		}
	}
	for (auto it : delay)
		pairDelay[it.first][host] = it.second;
	for (auto it : txDelay)
		pairTxDelay[it.first][host] = it.second;
	for (auto it : bw)
		pairBw[it.first->GetId()][host->GetId()] = it.second;
}

void CalculateRoutes(NodeContainer &n) {
	for (int i = 0; i < (int)n.GetN(); i++) {
		Ptr<Node> node = n.Get(i);
		if (node->GetNodeType() == 0)
			CalculateRoute(node);
	}
}

void SetRoutingEntries() {	//设置路由表条目
	// For each node.
	// Mapping destination to next hop for each node: <node, <dest, <nexthop0, ...> > >
	for (auto i = nextHop.begin(); i != nextHop.end(); i++) {
		Ptr<Node> node = i->first;
		auto &table = i->second;
		for (auto j = table.begin(); j != table.end(); j++) {
			// The destination node.
			Ptr<Node> dst = j->first;
			// The IP address of the dst.
			Ipv4Address dstAddr = dst->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal();
			// The next hops towards the dst.
			vector<Ptr<Node> > nexts = j->second;
			for (int k = 0; k < (int)nexts.size(); k++) {
				Ptr<Node> next = nexts[k];
				uint32_t interface = nbr2if[node][next].idx;
				if (node->GetNodeType()){
					if (node->GetNodeType() == 1){	//Tor switch
						//cout << "Setting Tbale Entry: " << node->GetId() << " to " << dst->GetId() << " Interface is " << interface << '\n';
						DynamicCast<TORSwitchNode>(node)->AddTableEntry(dstAddr, interface);
					}
					else if (node->GetNodeType() == 2) {//spine & core switch
						//cout << "Setting Tbale Entry: " << node->GetId() << " to " << dst->GetId() << " Interface is " << interface << '\n';
						DynamicCast<SwitchNode>(node)->AddTableEntry(dstAddr, interface);
					}
					else{// optical switch
						//cout << "Setting Tbale Entry: " << node->GetId() << " to " << dst->GetId() << " Interface is " << interface << '\n';
						DynamicCast<OSwitchNode>(node)->AddTableEntry(dstAddr, interface);
					}
				}
					
				else {
					//cout << "Setting Tbale Entry: " << node->GetId() << " to " << dst->GetId() << " Interface is " << interface << '\n';
					node->GetObject<RdmaDriver>()->m_rdma->AddTableEntry(dstAddr, interface);
				}
			}
		}
	}
}

void SetOpticalRoutes(){
	//Ptr<Node> node = OpticalSwitch;
	for (auto i = Tor2Server.begin(); i != Tor2Server.end(); i++){
		for (auto j = i->second.begin(); j!= i->second.end(); j++){
			Ptr<Node> dst = *j;
			Ipv4Address dstAddr = dst->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal();
			Ptr<Node> next = i->first;
			uint32_t interface = nbr2if[OpticalSwitch][next].idx;
			OpticalSwitch->AddTableEntry(dstAddr, interface);
		}
	}
}

// take down the link between a and b, and redo the routing
void TakeDownLink(NodeContainer n, Ptr<Node> a, Ptr<Node> b) {
	if (!nbr2if[a][b].up)
		return;
	// take down link between a and b
	nbr2if[a][b].up = nbr2if[b][a].up = false;
	nextHop.clear();
	CalculateRoutes(n);
	// clear routing tables
	for (uint32_t i = 0; i < n.GetN(); i++) {
		if (n.Get(i)->GetNodeType() == 1)
			DynamicCast<SwitchNode>(n.Get(i))->ClearTable();
		else
			n.Get(i)->GetObject<RdmaDriver>()->m_rdma->ClearTable();
	}
	DynamicCast<QbbNetDevice>(a->GetDevice(nbr2if[a][b].idx))->TakeDown();
	DynamicCast<QbbNetDevice>(b->GetDevice(nbr2if[b][a].idx))->TakeDown();
	// reset routing table
	SetRoutingEntries();

	// redistribute qp on each host
	for (uint32_t i = 0; i < n.GetN(); i++) {
		if (n.Get(i)->GetNodeType() == 0)
			n.Get(i)->GetObject<RdmaDriver>()->m_rdma->RedistributeQp();
	}
}

uint64_t get_nic_rate(NodeContainer &n) {
	for (uint32_t i = 0; i < n.GetN(); i++)
		if (n.Get(i)->GetNodeType() == 0)
			return DynamicCast<QbbNetDevice>(n.Get(i)->GetDevice(1))->GetDataRate().GetBitRate();
	std::cout << "get nic rate error!\n";
	exit(-1);
}

//获取ToR交换机的吞吐率情况
void PrintResults(std::map<uint32_t, NetDeviceContainer> ToR, uint32_t numToRs, double delay) {
	//获取ToR交换机的带宽结果
	for (uint32_t i = 0; i < numToRs; i++) {
		double throughputTotal = 0;
		uint64_t torBuffer = 0;
		double power;
		for (uint32_t j = 0; j < ToR[i].GetN(); j++) {
			Ptr<QbbNetDevice> nd = DynamicCast<QbbNetDevice>(ToR[i].Get(j));
//			uint64_t txBytes = nd->getTxBytes();
			uint64_t txBytes = nd->GetQueue()->getTxBytes();
			double rxBytes = nd->getNumRxBytes();

			uint64_t qlen = nd->GetQueue()->GetNBytesTotal();
			uint64_t bw = nd->GetDataRate().GetBitRate(); //maxRtt

			torBuffer += qlen;
			double throughput = double(txBytes * 8) / delay;
			if (j == 2) { //  ToDo. very ugly hardcode here specific to the burst evaluation scenario where 4 is the receiver in flow-burstExp.txt.
				throughputTotal += throughput;
				power = (rxBytes * 8.0 / delay) * (qlen + bw * maxRtt * 1e-9) / (bw * (bw * maxRtt * 1e-9));
			}
			std::cout << "ToR " << i << " Port " << j << " throughput " << throughput << " txBytes " << txBytes << " qlen " << qlen << " time " << Simulator::Now().GetSeconds() << " normpower " << power << std::endl;
		}
		std::cout << "ToR " << i << " Total " << 0 << " throughput " << throughputTotal << " buffer " << torBuffer <<  " time " << Simulator::Now().GetSeconds() << std::endl;
	}
	Simulator::Schedule(Seconds(delay), PrintResults, ToR, numToRs, delay);
}

//获取每一个服务器Src的吞吐率？
void PrintResultsFlow(std::map<uint32_t, NetDeviceContainer> Src, uint32_t numFlows, double delay) {
	//获取流的完成结果
	for (uint32_t i = 0; i < numFlows; i++) {
		double throughputTotal = 0;

		for (uint32_t j = 0; j < Src[i].GetN(); j++) {
			Ptr<QbbNetDevice> nd = DynamicCast<QbbNetDevice>(Src[i].Get(j));
//			uint64_t txBytes = nd->getTxBytes();
			uint64_t txBytes = nd->getNumTxBytes();

			uint64_t qlen = nd->GetQueue()->GetNBytesTotal();
			double throughput = double(txBytes * 8) / delay;
			throughputTotal += throughput;
			// std::cout << "Src " << i << " Port " << j << " throughput "<< throughput << " txBytes " << txBytes << " qlen " << qlen << " time " << Simulator::Now().GetSeconds() << std::endl;
		}
		std::cout << "Src " << i << " Total " << 0 << " throughput " << throughputTotal <<  " time " << Simulator::Now().GetSeconds() << std::endl;
	}
	Simulator::Schedule(Seconds(delay), PrintResultsFlow, Src, numFlows, delay);
}

void show_nexthops(){
	for (const auto& outerPair : nextHop) {
        Ptr<Node> key1 = outerPair.first;
        const auto& innerMap = outerPair.second;
        
        std::cout << "Key1 (Node ID): " << key1->GetId() << std::endl;
        
        for (const auto& innerPair : innerMap) {
            Ptr<Node> key2 = innerPair.first;
            const auto& nodeVector = innerPair.second;
            
            std::cout << "  Key2 (Node ID): " << key2->GetId() << std::endl;
            std::cout << "    Vector of Nodes: ";
            
            for (const auto& nodePtr : nodeVector) {
                std::cout << nodePtr->GetId() << " ";
            }
            std::cout << std::endl;
        }
    }
}

void PrintNbr2If() {
    for (const auto& outerPair : nbr2if) {
        Ptr<Node> key1 = outerPair.first;
        const auto& innerMap = outerPair.second;
        
        std::cout << "Key1 (Node ID): " << key1->GetId() << std::endl;
        
        for (const auto& innerPair : innerMap) {
            Ptr<Node> key2 = innerPair.first;
            const Interface& iface = innerPair.second;
            
            std::cout << "  Key2 (Node ID): " << key2->GetId() << std::endl;
            std::cout << "    Interface Info: " << iface.idx << std::endl;
        }
    }
}

int main(int argc, char* argv[])
{
    clock_t begint, endt; // 记录时间，clock_t 为某种 long 类型
    begint = clock(); // clock() 函数获取当前时间
	std::ifstream conf; // 声明一个文件输入流 conf，用于读取配置文件

	// 定义两个 bool 类型变量，用于启用不同的 TCP 拥塞控制算法
	bool wien = true; // wien enables PowerTCP. 
	bool delayWien = false; // delayWien enables Theta-PowerTCP (delaypowertcp) 

    uint32_t algorithm = 3; // 拥塞控制算法类型
	uint32_t windowCheck = 1; // 窗口检查策略
	std::string confFile = "/home/ss/CRANE-ns3/simulator/ns-3.39/examples/PowerTCP/config-first_ss.txt"; // 定义配置文件的路径
	std::cout << confFile; // 输出检查，用于调试

	CommandLine cmd; // 创建一个 CommandLine 对象，用于解析命令行参数
	// 将命令行参数与程序中的变量绑定，使得可以通过命令行来覆盖默认配置
	cmd.AddValue("conf", "config file path", confFile);
	cmd.AddValue("wien", "enable wien --> wien enables PowerTCP.", wien);
	cmd.AddValue("delayWien", "enable wien delay --> delayWien enables Theta-PowerTCP (delaypowertcp) ", delayWien);
	cmd.AddValue ("algorithm", "specify CC mode. This is added for my convinience. I prefer cmd rather than parsing files.", algorithm);
	cmd.AddValue("windowCheck", "windowCheck", windowCheck);
	cmd.Parse (argc, argv); // 解析命令行参数，将其存储到之前定义的变量中

	conf.open(confFile.c_str()); // 打开指定的配置文件
    while (!conf.eof()) // 开始循环，直到文件末尾
	{
		// 声明一个字符串 key ，并从配置文件中读取一行
		std::string key;
		conf >> key; // 读取下一个字符串

		// 是否开启 QCN
		if (key.compare("ENABLE_QCN") == 0) // 字符串是否匹配，== 0 则表示匹配
		{
			uint32_t v;
			conf >> v;
			enable_qcn = v; // 将 QCN 开启与否的情况赋值给 enable_qcn
			if (enable_qcn)
				std::cout << "ENABLE_QCN\t\t\t" << "Yes" << "\n"; // 打印输出，便于检查
			else
				std::cout << "ENABLE_QCN\t\t\t" << "No" << "\n";
		}
		// 是否限制目标速率
		else if (key.compare("CLAMP_TARGET_RATE") == 0) 
		{
			uint32_t v;
			conf >> v;
			clamp_target_rate = v;
			if (clamp_target_rate)
				std::cout << "CLAMP_TARGET_RATE\t\t" << "Yes" << "\n";
			else
				std::cout << "CLAMP_TARGET_RATE\t\t" << "No" << "\n";
		}
		else if (key.compare("PAUSE_TIME") == 0)
		{
			double v;
			conf >> v;
			pause_time = v;
			std::cout << "PAUSE_TIME\t\t\t" << pause_time << "\n";
		}
		else if (key.compare("DATA_RATE") == 0)
		{
			std::string v;
			conf >> v;
			data_rate = v;
			std::cout << "DATA_RATE\t\t\t" << data_rate << "\n";
		}
		else if (key.compare("LINK_DELAY") == 0)
		{
			std::string v;
			conf >> v;
			link_delay = v;
			std::cout << "LINK_DELAY\t\t\t" << link_delay << "\n";
		}
		else if (key.compare("PACKET_PAYLOAD_SIZE") == 0) // 数据包负载大小
		{
			uint32_t v;
			conf >> v;
			packet_payload_size = v;
			std::cout << "PACKET_PAYLOAD_SIZE\t\t" << packet_payload_size << "\n";
		}
		else if (key.compare("L2_CHUNK_SIZE") == 0) // 链路层数据块大小？
		{
			uint32_t v;
			conf >> v;
			l2_chunk_size = v;
			std::cout << "L2_CHUNK_SIZE\t\t\t" << l2_chunk_size << "\n";
		}
		else if (key.compare("L2_ACK_INTERVAL") == 0) //数据层 ACK 发送间隔
		{
			uint32_t v;
			conf >> v;
			l2_ack_interval = v;
			std::cout << "L2_ACK_INTERVAL\t\t\t" << l2_ack_interval << "\n";
		}
		else if (key.compare("L2_BACK_TO_ZERO") == 0)
		{
			uint32_t v;
			conf >> v;
			l2_back_to_zero = v;
			if (l2_back_to_zero)
				std::cout << "L2_BACK_TO_ZERO\t\t\t" << "Yes" << "\n";
			else
				std::cout << "L2_BACK_TO_ZERO\t\t\t" << "No" << "\n";
		}
		else if (key.compare("TOPOLOGY_FILE") == 0)
		{
			std::string v;
			conf >> v;
			topology_file = v;
			std::cout << "TOPOLOGY_FILE\t\t\t" << topology_file << "\n";
		}
		else if (key.compare("FLOW_FILE") == 0)
		{
			std::string v;
			conf >> v;
			flow_file = v;
			std::cout << "FLOW_FILE\t\t\t" << flow_file << "\n";
		}
		else if (key.compare("TRACE_FILE") == 0)
		{
			std::string v;
			conf >> v;
			trace_file = v;
			std::cout << "TRACE_FILE\t\t\t" << trace_file << "\n";
		}
		else if (key.compare("TRACE_OUTPUT_FILE") == 0)
		{
			std::string v;
			conf >> v;
			trace_output_file = v;
			if (argc > 2)
			{
				trace_output_file = trace_output_file + std::string(argv[2]);
			}
			std::cout << "TRACE_OUTPUT_FILE\t\t" << trace_output_file << "\n";
		}
		else if (key.compare("SIMULATOR_STOP_TIME") == 0)
		{
			double v;
			conf >> v;
			simulator_stop_time = v;
			std::cout << "SIMULATOR_STOP_TIME\t\t" << simulator_stop_time << "\n";
		}
		else if (key.compare("ALPHA_RESUME_INTERVAL") == 0)
		{
			double v;
			conf >> v;
			alpha_resume_interval = v;
			std::cout << "ALPHA_RESUME_INTERVAL\t\t" << alpha_resume_interval << "\n";
		}
		else if (key.compare("RP_TIMER") == 0)
		{
			double v;
			conf >> v;
			rp_timer = v;
			std::cout << "RP_TIMER\t\t\t" << rp_timer << "\n";
		}
		else if (key.compare("EWMA_GAIN") == 0)
		{
			double v;
			conf >> v;
			ewma_gain = v;
			std::cout << "EWMA_GAIN\t\t\t" << ewma_gain << "\n";
		}
		else if (key.compare("FAST_RECOVERY_TIMES") == 0)
		{
			uint32_t v;
			conf >> v;
			fast_recovery_times = v;
			std::cout << "FAST_RECOVERY_TIMES\t\t" << fast_recovery_times << "\n";
		}
		else if (key.compare("RATE_AI") == 0)
		{
			std::string v;
			conf >> v;
			rate_ai = v;
			std::cout << "RATE_AI\t\t\t\t" << rate_ai << "\n";
		}
		else if (key.compare("RATE_HAI") == 0)
		{
			std::string v;
			conf >> v;
			rate_hai = v;
			std::cout << "RATE_HAI\t\t\t" << rate_hai << "\n";
		}
		else if (key.compare("ERROR_RATE_PER_LINK") == 0)
		{
			double v;
			conf >> v;
			error_rate_per_link = v;
			std::cout << "ERROR_RATE_PER_LINK\t\t" << error_rate_per_link << "\n";
		}
		else if (key.compare("CC_MODE") == 0) {
			conf >> cc_mode;
			std::cout << "CC_MODE\t\t" << cc_mode << '\n';
		} else if (key.compare("RATE_DECREASE_INTERVAL") == 0) {
			double v;
			conf >> v;
			rate_decrease_interval = v;
			std::cout << "RATE_DECREASE_INTERVAL\t\t" << rate_decrease_interval << "\n";
		} else if (key.compare("MIN_RATE") == 0) {
			conf >> min_rate;
			std::cout << "MIN_RATE\t\t" << min_rate << "\n";
		} else if (key.compare("FCT_OUTPUT_FILE") == 0) {
			conf >> fct_output_file;
			std::cout << "FCT_OUTPUT_FILE\t\t" << fct_output_file << '\n';
		} else if (key.compare("HAS_WIN") == 0) {
			conf >> has_win;
			std::cout << "HAS_WIN\t\t" << has_win << "\n";
		} else if (key.compare("GLOBAL_T") == 0) {
			conf >> global_t;
			std::cout << "GLOBAL_T\t\t" << global_t << '\n';
		} else if (key.compare("MI_THRESH") == 0) {
			conf >> mi_thresh;
			std::cout << "MI_THRESH\t\t" << mi_thresh << '\n';
		} else if (key.compare("VAR_WIN") == 0) {
			uint32_t v;
			conf >> v;
			var_win = v;
			std::cout << "VAR_WIN\t\t" << v << '\n';
		} else if (key.compare("FAST_REACT") == 0) {
			uint32_t v;
			conf >> v;
			fast_react = v;
			std::cout << "FAST_REACT\t\t" << v << '\n';
		} else if (key.compare("U_TARGET") == 0) {
			conf >> u_target;
			std::cout << "U_TARGET\t\t" << u_target << '\n';
		} else if (key.compare("INT_MULTI") == 0) {
			conf >> int_multi;
			std::cout << "INT_MULTI\t\t\t\t" << int_multi << '\n';
		} else if (key.compare("RATE_BOUND") == 0) {
			uint32_t v;
			conf >> v;
			rate_bound = v;
			std::cout << "RATE_BOUND\t\t" << rate_bound << '\n';
		} else if (key.compare("ACK_HIGH_PRIO") == 0) {
			conf >> ack_high_prio;
			std::cout << "ACK_HIGH_PRIO\t\t" << ack_high_prio << '\n';
		} else if (key.compare("DCTCP_RATE_AI") == 0) {
			conf >> dctcp_rate_ai;
			std::cout << "DCTCP_RATE_AI\t\t\t\t" << dctcp_rate_ai << "\n";
		} else if (key.compare("PFC_OUTPUT_FILE") == 0) {
			conf >> pfc_output_file;
			std::cout << "PFC_OUTPUT_FILE\t\t\t\t" << pfc_output_file << '\n';
		} else if (key.compare("LINK_DOWN") == 0) {
			conf >> link_down_time >> link_down_A >> link_down_B;
			std::cout << "LINK_DOWN\t\t\t\t" << link_down_time << ' ' << link_down_A << ' ' << link_down_B << '\n';
		} else if (key.compare("ENABLE_TRACE") == 0) {
			conf >> enable_trace;
			std::cout << "ENABLE_TRACE\t\t\t\t" << enable_trace << '\n';
		} else if (key.compare("KMAX_MAP") == 0) {
			int n_k ;
			conf >> n_k;
			std::cout << "KMAX_MAP\t\t\t\t";
			for (int i = 0; i < n_k; i++) {
				uint64_t rate;
				uint32_t k;
				conf >> rate >> k;
				rate2kmax[rate] = k;
				std::cout << ' ' << rate << ' ' << k;
			}
			std::cout << '\n';
		} else if (key.compare("KMIN_MAP") == 0) {
			int n_k ;
			conf >> n_k;
			std::cout << "KMIN_MAP\t\t\t\t";
			for (int i = 0; i < n_k; i++) {
				uint64_t rate;
				uint32_t k;
				conf >> rate >> k;
				rate2kmin[rate] = k;
				std::cout << ' ' << rate << ' ' << k;
			}
			std::cout << '\n';
		} else if (key.compare("PMAX_MAP") == 0) {
			int n_k ;
			conf >> n_k;
			std::cout << "PMAX_MAP\t\t\t\t";
			for (int i = 0; i < n_k; i++) {
				uint64_t rate;
				double p;
				conf >> rate >> p;
				rate2pmax[rate] = p;
				std::cout << ' ' << rate << ' ' << p;
			}
			std::cout << '\n';
		} else if (key.compare("BUFFER_SIZE") == 0) {
			conf >> buffer_size;
			std::cout << "BUFFER_SIZE\t\t\t\t" << buffer_size << '\n';
		} else if (key.compare("QLEN_MON_FILE") == 0) {
			conf >> qlen_mon_file;
			std::cout << "QLEN_MON_FILE\t\t\t\t" << qlen_mon_file << '\n';
		} else if (key.compare("QLEN_MON_START") == 0) {
			conf >> qlen_mon_start;
			std::cout << "QLEN_MON_START\t\t\t\t" << qlen_mon_start << '\n';
		} else if (key.compare("QLEN_MON_END") == 0) {
			conf >> qlen_mon_end;
			std::cout << "QLEN_MON_END\t\t\t\t" << qlen_mon_end << '\n';
		} else if (key.compare("MULTI_RATE") == 0) {
			int v;
			conf >> v;
			multi_rate = v;
			std::cout << "MULTI_RATE\t\t\t\t" << multi_rate << '\n';
		} else if (key.compare("SAMPLE_FEEDBACK") == 0) {
			int v;
			conf >> v;
			sample_feedback = v;
			std::cout << "SAMPLE_FEEDBACK\t\t\t\t" << sample_feedback << '\n';
		} else if (key.compare("PINT_LOG_BASE") == 0) {
			conf >> pint_log_base;
			std::cout << "PINT_LOG_BASE\t\t\t\t" << pint_log_base << '\n';
		} else if (key.compare("PINT_PROB") == 0) {
			conf >> pint_prob;
			std::cout << "PINT_PROB\t\t\t\t" << pint_prob << '\n';
		} else if (key.compare("OSWITCH_CONFIGURE_DELAY") == 0) {
			uint64_t v;
			conf >> v;
			optical_reconfigure_delay = v;
			std::cout << "OSWITCH_CONFIGURE_DELAY\t\t\t\t" << optical_reconfigure_delay << '\n';
		}else if (key.compare("OSWITCH_PACKET_RATE") == 0) {
			uint64_t v;
			conf >> v;
			optical_pswitch_rate = v;
			std::cout << "OSWITCH_PACKET_RATE\t\t\t\t" << optical_pswitch_rate << '\n';
		}else if (key.compare("OSIWTCH_OPTICAL_RATE") == 0) {
			uint64_t v;
			conf >> v;
			optical_oswitch_rate = v;
			std::cout << "OSIWTCH_OPTICAL_RATE\t\t\t\t" << optical_oswitch_rate << '\n';
		}else if (key.compare("OSWITCH_DURATION_EXTENTIONRATE") == 0) {
			double v;
			conf >> v;
			optical_theta = v;
			std::cout << "OSWITCH_DURATION_EXTENTIONRATE\t\t\t\t" << optical_theta << '\n';
		}else if (key.compare("OPTICAL_FILE") == 0) {
			conf >> optical_conf_file;
			std:: cout << "OPTICAL_FILE\t\t\t\t" << optical_conf_file << '\n';
		}
		fflush(stdout); // 清空标准输出缓冲区
	}
    conf.close(); // 结束循环，关闭配置文件

	// 将命令行参数或配置文件中的值赋给全局变量，覆盖默认值
    // cc_mode = algorithm; // overrides configuration file
	cc_mode = 0;
	has_win = windowCheck; // overrides configuration file
	var_win = windowCheck; // overrides configuration file

	// 设置 NS-3 默认配置，应用全局参数
    Config::SetDefault("ns3::QbbNetDevice::PauseTime", UintegerValue(pause_time));
	Config::SetDefault("ns3::QbbNetDevice::QcnEnabled", BooleanValue(enable_qcn));

    // 设置 IntHop 的多路径参数，影响拥塞控制算法
	IntHop::multi = int_multi;

	// 根据拥塞控制模式设置头部类型
	if (cc_mode == 7) // timely, use ts
		// IntHeader::mode = IntHeader::TS;
		IntHeader::mode = IntHeader::NONE;
	else if (cc_mode == 3) // hpcc, powertcp, use int
		// IntHeader::mode = IntHeader::NORMAL;
		IntHeader::mode = IntHeader::NONE;
	else if (cc_mode == 10) // hpcc-pint
		// IntHeader::mode = IntHeader::PINT;
		IntHeader::mode = IntHeader::NONE;
	else // others, no extra header
		IntHeader::mode = IntHeader::NONE;

	// 如果使用 Pint，设置其日志基数和字节数
	if (cc_mode == 10) {
		// Pint::set_log_base(pint_log_base);
		// IntHeader::pint_bytes = Pint::get_n_bytes();
	}

    topof.open(topology_file.c_str()); // 打开拓扑文件
	flowf.open(flow_file.c_str()); // 打开流量文件
	// 定义节点数量、OCS 数量、EPS 数量、ToR 交换机数量和链路数量
	// 其中，EPS 数量包含了 ToR 交换机数量 
	uint32_t node_num, oswitch_num, switch_num, tors, link_num; // 139 1 10 4 148

	// 输入拓扑情况
    topof >> node_num >> oswitch_num >>  switch_num >> tors >> link_num; 
    std::cout << node_num << " " << switch_num << " " << tors <<  " " << link_num << std::endl;
	// 输入流量情况
    flowf >> flow_num;

    NodeContainer serverNodes; // 服务器节点
	NodeContainer torNodes; // ToR 交换机节点
	NodeContainer spineNodes; // 
	NodeContainer switchNodes; // EPS 节点
    NodeContainer oswitchNodes; // OCS 节点
	NodeContainer allNodes; // 

    std::vector<uint32_t> node_type(node_num, 0); // 创建一个向数组，初始化节点类型为 0
    std::cout << "switch_num " << switch_num << std::endl;

	/*
		假设输入为：
		139 1 10 4 148
		129 130 131 132 133 134 135 136 137 138
		则有 10 台 EPS，且 10 台中的前 4 台为 ToR 交换机
	*/
    for (uint32_t i = 0; i < switch_num; i++) {
		uint32_t sid;
		topof >> sid;
		std::cout << "sid " << sid << std::endl;
		switchNumToId[i] = sid; // map 类型，键从 0 开始编号，值为 EPS 的id
		switchIdToNum[sid] = i; // map 类型，键为 EPS 的id，值为 num 编号
		if (i < tors) {
			node_type[sid] = 1; // 属于 ToR 交换机
		}
		else
			node_type[sid] = 2; // 属于 spine 或 core 交换机
	}

    node_type[0] = 3; // for optical switch. Optical switch must be front of Tor.
    //Ptr<OSwitchNode> ocs = CreateObject<OSwitchNode>();

    for (uint32_t i = 0; i < node_num; i++) { // 遍历所有节点
		if (node_type[i] == 0) { // 服务器节点
			Ptr<Node> node = CreateObject<Node>();
			n.Add(node);
			allNodes.Add(node);
			serverNodes.Add(node);
			//cout << "Created server, node id is: " << node->GetId() << '\n';
		}
        else if (node_type[i] == 3){ // OCS 节点
            Ptr<OSwitchNode> ocs = CreateObject<OSwitchNode>();
            OpticalSwitch = ocs;
            n.Add(ocs);
			allNodes.Add(ocs);
            oswitchNodes.Add(ocs);
            ocs->SetNodeType(3);
			//cout << "Created Oswitch, node id is: " << ocs->GetId() << '\n';
        }
        else if (node_type[i] == 1){ // ToR 交换机节点
            Ptr<TORSwitchNode> tsw = CreateObject<TORSwitchNode>();
            n.Add(tsw);
            switchNodes.Add(tsw);
            allNodes.Add(tsw);
            tsw->SetAttribute("EcnEnabled", BooleanValue(enable_qcn));
            torNodes.Add(tsw);
            tsw->SetNodeType(1);
			tsw->SetOCS(OpticalSwitch);
			// tsw->Add_m_server_num();
			// tsw->Add_m_server_num();
			//cout << "Created ToR switch, node id is: " << tsw->GetId() << '\n';
        }
		else { // spine 和 core 交换机节点
			Ptr<SwitchNode> sw = CreateObject<SwitchNode>();
			n.Add(sw);
			switchNodes.Add(sw);
			allNodes.Add(sw);
			sw->SetAttribute("EcnEnabled", BooleanValue(enable_qcn));
			spineNodes.Add(sw);
			sw->SetNodeType(2);
			//cout << "Created switch, node id is: " << sw->GetId() << '\n';
		}
	}

	cout << "Created " << n.GetN() << " Nodes\n"; // 输出创建的节点数量

    NS_LOG_INFO("Create nodes."); // 记录日志信息，标记节点创建完成

	// 创建互联网堆栈和全局路由助手，并将其安装到所有节点上
    InternetStackHelper internet;
	Ipv4GlobalRoutingHelper globalRoutingHelper;
	internet.SetRoutingHelper (globalRoutingHelper);
	//cout << "Preaparing to inetnert.Install(n).\n";
	internet.Install(n);
	//cout << "Preaparing to Assign IP to each server.\n";
    
	// 为每个服务器节点分配 IP 地址
    for (uint32_t i = 0; i < node_num; i++) {
		if (n.Get(i)->GetNodeType() == 0) { // is server
			serverAddress.resize(i + 1);
			serverAddress[i] = node_id_to_ip(i);
		}
	}

    NS_LOG_INFO("Create channels."); // 记录日志信息，标记通道创建开始

	// 创建并配置一个速率错误模型，用于模拟链路错误率
    Ptr<RateErrorModel> rem = CreateObject<RateErrorModel>();
	Ptr<UniformRandomVariable> uv = CreateObject<UniformRandomVariable>();
	rem->SetRandomVariable(uv);
	uv->SetStream(50);
	rem->SetAttribute("ErrorRate", DoubleValue(error_rate_per_link));
	rem->SetAttribute("ErrorUnit", StringValue("ERROR_UNIT_PACKET"));

    FILE *pfc_file = fopen(pfc_output_file.c_str(), "w");

    QbbHelper qbb;
	Ipv4AddressHelper ipv4;

	// 循环读取链路信息，包括源节点、目标节点、数据速率、链路延迟和错误率
    for (uint32_t i = 0; i < link_num; i++)
	{
		uint32_t src, dst;
		std::string data_rate, link_delay;
		double error_rate;
		topof >> src >> dst >> data_rate >> link_delay >> error_rate;
		// n.GetN函数，用于获取nodecontainer的size
		std::cout << src << " " << dst << " " << n.GetN() << " " << data_rate << " " << link_delay << " " << error_rate << std::endl;
		Ptr<Node> snode = n.Get(src), dnode = n.Get(dst); // 获取源节点和目标节点的指针

		qbb.SetDeviceAttribute("DataRate", StringValue(data_rate)); // 设置链路速率
		qbb.SetChannelAttribute("Delay", StringValue(link_delay)); // 设置链路延迟
		if (error_rate > 0)
		{
			Ptr<RateErrorModel> rem = CreateObject<RateErrorModel>();
			Ptr<UniformRandomVariable> uv = CreateObject<UniformRandomVariable>();
			rem->SetRandomVariable(uv);
			uv->SetStream(50);
			rem->SetAttribute("ErrorRate", DoubleValue(error_rate));
			rem->SetAttribute("ErrorUnit", StringValue("ERROR_UNIT_PACKET"));
			qbb.SetDeviceAttribute("ReceiveErrorModel", PointerValue(rem));
		}
		else
		{
			qbb.SetDeviceAttribute("ReceiveErrorModel", PointerValue(rem));
		}

		fflush(stdout); // 刷新标准输出缓冲区

		// Assigne server IP
		// Note: this should be before the automatic assignment below (ipv4.Assign(d)),
		// because we want our IP to be the primary IP (first in the IP address list),
		// so that the global routing is based on our IP
		NetDeviceContainer d = qbb.Install(snode, dnode);// Helper 帮助两个 Node 构建netdevice，并且连接到channel中

		// 如果源节点是服务器，则为其分配 IP 地址
		if (snode->GetNodeType() == 0) { 
			Ptr<Ipv4> ipv4 = snode->GetObject<Ipv4>();
			ipv4->AddInterface(d.Get(0));
			ipv4->AddAddress(1, Ipv4InterfaceAddress(serverAddress[src], Ipv4Mask(0xff000000)));
		}
		// 如果目标节点是服务器，同样为其分配 IP 地址
		if (dnode->GetNodeType() == 0) { 
			Ptr<Ipv4> ipv4 = dnode->GetObject<Ipv4>();
			ipv4->AddInterface(d.Get(1));
			ipv4->AddAddress(1, Ipv4InterfaceAddress(serverAddress[dst], Ipv4Mask(0xff000000)));
		}

		if (!snode->GetNodeType()) {	//src为server
			sourceNodes[src].Add(DynamicCast<QbbNetDevice>(d.Get(0)));
		}

		if (!snode->GetNodeType() && dnode->GetNodeType()) {	//从服务器到交换机,拓扑文件中不存在交换机->服务器的拓扑，只考虑这一种情况
			switchDown[switchIdToNum[dst]].Add(DynamicCast<QbbNetDevice>(d.Get(1)));
			Tor2Server[dnode].push_back(snode);

			if(dnode->GetNodeType() == 1){//更新ToR的server number
				Ptr<TORSwitchNode> tsw = DynamicCast<TORSwitchNode>(dnode);
				tsw->Add_m_server_num();
			}
		}


		if (snode->GetNodeType() && dnode->GetNodeType()) {		//交换机到交换机
			switchToSwitchInterfaces.Add(d);
			switchUp[switchIdToNum[src]].Add(DynamicCast<QbbNetDevice>(d.Get(0)));
			switchUp[switchIdToNum[dst]].Add(DynamicCast<QbbNetDevice>(d.Get(1)));
			switchToSwitch[src][dst].push_back(DynamicCast<QbbNetDevice>(d.Get(0)));
			switchToSwitch[src][dst].push_back(DynamicCast<QbbNetDevice>(d.Get(1)));
		}

		//cout << "Going to create a graph of the topology, setting nbr2if.\n";

		// used to create a graph of the topology
		nbr2if[snode][dnode].idx = DynamicCast<QbbNetDevice>(d.Get(0))->GetIfIndex();
		nbr2if[snode][dnode].up = true;
		nbr2if[snode][dnode].delay = DynamicCast<QbbChannel>(DynamicCast<QbbNetDevice>(d.Get(0))->GetChannel())->GetDelay().GetTimeStep();
		nbr2if[snode][dnode].bw = DynamicCast<QbbNetDevice>(d.Get(0))->GetDataRate().GetBitRate();
		nbr2if[dnode][snode].idx = DynamicCast<QbbNetDevice>(d.Get(1))->GetIfIndex();
		nbr2if[dnode][snode].up = true;
		nbr2if[dnode][snode].delay = DynamicCast<QbbChannel>(DynamicCast<QbbNetDevice>(d.Get(1))->GetChannel())->GetDelay().GetTimeStep();
		nbr2if[dnode][snode].bw = DynamicCast<QbbNetDevice>(d.Get(1))->GetDataRate().GetBitRate();

		//cout << "nbr2if set.\n";

		// This is just to set up the connectivity between nodes. The IP addresses are useless
		// char ipstring[16];
		std::stringstream ipstring;
		ipstring << "10." << i / 254 + 1 << "." << i % 254 + 1 << ".0";
		// sprintf(ipstring, "10.%d.%d.0", i / 254 + 1, i % 254 + 1);
		ipv4.SetBase(ipstring.str().c_str(), "255.255.255.0");
		ipv4.Assign(d);

		// setup PFC trace
		// DynamicCast<QbbNetDevice>(d.Get(0))->TraceConnectWithoutContext("QbbPfc", MakeBoundCallback (&get_pfc, pfc_file, DynamicCast<QbbNetDevice>(d.Get(0))));
		// DynamicCast<QbbNetDevice>(d.Get(1))->TraceConnectWithoutContext("QbbPfc", MakeBoundCallback (&get_pfc, pfc_file, DynamicCast<QbbNetDevice>(d.Get(1))));
	}

    nic_rate = get_nic_rate(n);
	// cout << "get_nic_rate done.\n";
    // The switch mmu runs Dynamic Thresholds (DT) by default.
	//cout << "The switch mmu begin.\n";
    for (uint32_t i = 0; i < node_num; i++) {
		if (n.Get(i)->GetNodeType() && n.Get(i)->GetNodeType() !=3) { // is switch, and not optical switch
            if(n.Get(i)->GetNodeType() == 1) { // tor switch
                Ptr<TORSwitchNode> sw = DynamicCast<TORSwitchNode>(n.Get(i));
                uint32_t shift = 3; // by default 1/8
                double alpha = 1.0 / 8;
                sw->m_mmu->SetAlphaIngress(alpha);
                sw->m_mmu->SetAlphaEgress(UINT16_MAX);
                uint64_t totalHeadroom = 0;
                for (uint32_t j = 1; j < sw->GetNDevices(); j++) {

                    for (uint32_t qu = 0; qu < 8; qu++) {
                        Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(sw->GetDevice(j));
                        // set ecn
                        uint64_t rate = dev->GetDataRate().GetBitRate();
                        NS_ASSERT_MSG(rate2kmin.find(rate) != rate2kmin.end(), "must set kmin for each link speed");
                        NS_ASSERT_MSG(rate2kmax.find(rate) != rate2kmax.end(), "must set kmax for each link speed");
                        NS_ASSERT_MSG(rate2pmax.find(rate) != rate2pmax.end(), "must set pmax for each link speed");
                        sw->m_mmu->ConfigEcn(j, rate2kmin[rate], rate2kmax[rate], rate2pmax[rate]);
                        // set pfc
                        uint64_t delay = DynamicCast<QbbChannel>(dev->GetChannel())->GetDelay().GetTimeStep();
                        uint32_t headroom = rate * delay / 8 / 1000000000 * 3;

                        sw->m_mmu->SetHeadroom(headroom, j, qu);
                        totalHeadroom += headroom;
                    }

                }
                sw->m_mmu->SetBufferPool(buffer_size * 1024 * 1024);
                sw->m_mmu->SetIngressPool(buffer_size * 1024 * 1024 - totalHeadroom);
                sw->m_mmu->SetEgressLosslessPool(buffer_size * 1024 * 1024);
                sw->m_mmu->node_id = sw->GetId();
            }
            else{
                Ptr<SwitchNode> sw = DynamicCast<SwitchNode>(n.Get(i));
                uint32_t shift = 3; // by default 1/8
                double alpha = 1.0 / 8;
                sw->m_mmu->SetAlphaIngress(alpha);
                sw->m_mmu->SetAlphaEgress(UINT16_MAX);
                uint64_t totalHeadroom = 0;
                for (uint32_t j = 1; j < sw->GetNDevices(); j++) {

                    for (uint32_t qu = 0; qu < 8; qu++) {
                        Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(sw->GetDevice(j));
                        // set ecn
                        uint64_t rate = dev->GetDataRate().GetBitRate();
                        NS_ASSERT_MSG(rate2kmin.find(rate) != rate2kmin.end(), "must set kmin for each link speed");
                        NS_ASSERT_MSG(rate2kmax.find(rate) != rate2kmax.end(), "must set kmax for each link speed");
                        NS_ASSERT_MSG(rate2pmax.find(rate) != rate2pmax.end(), "must set pmax for each link speed");
                        sw->m_mmu->ConfigEcn(j, rate2kmin[rate], rate2kmax[rate], rate2pmax[rate]);
                        // set pfc
                        uint64_t delay = DynamicCast<QbbChannel>(dev->GetChannel())->GetDelay().GetTimeStep();
                        uint32_t headroom = rate * delay / 8 / 1000000000 * 3;

                        sw->m_mmu->SetHeadroom(headroom, j, qu);
                        totalHeadroom += headroom;
                    }

                }
                sw->m_mmu->SetBufferPool(buffer_size * 1024 * 1024);
                sw->m_mmu->SetIngressPool(buffer_size * 1024 * 1024 - totalHeadroom);
                sw->m_mmu->SetEgressLosslessPool(buffer_size * 1024 * 1024);
                sw->m_mmu->node_id = sw->GetId();
            }
		}
	}
	//cout << "The switch mmu done.\n";
#if ENABLE_QP
	FILE *fct_output = fopen(fct_output_file.c_str(), "w");
	//
	// install RDMA driver 	
	//
	for (uint32_t i = 0; i < node_num; i++) {
		if (n.Get(i)->GetNodeType() == 0) { // is server
			// create RdmaHw
			Ptr<RdmaHw> rdmaHw = CreateObject<RdmaHw>();
			rdmaHw->SetAttribute("ClampTargetRate", BooleanValue(clamp_target_rate));
			rdmaHw->SetAttribute("AlphaResumInterval", DoubleValue(alpha_resume_interval));
			rdmaHw->SetAttribute("RPTimer", DoubleValue(rp_timer));
			rdmaHw->SetAttribute("FastRecoveryTimes", UintegerValue(fast_recovery_times));
			rdmaHw->SetAttribute("EwmaGain", DoubleValue(ewma_gain));
			rdmaHw->SetAttribute("RateAI", DataRateValue(DataRate(rate_ai)));
			rdmaHw->SetAttribute("RateHAI", DataRateValue(DataRate(rate_hai)));
			rdmaHw->SetAttribute("L2BackToZero", BooleanValue(l2_back_to_zero));
			rdmaHw->SetAttribute("L2ChunkSize", UintegerValue(l2_chunk_size));
			rdmaHw->SetAttribute("L2AckInterval", UintegerValue(l2_ack_interval));
			rdmaHw->SetAttribute("CcMode", UintegerValue(cc_mode));
			rdmaHw->SetAttribute("RateDecreaseInterval", DoubleValue(rate_decrease_interval));
			rdmaHw->SetAttribute("MinRate", DataRateValue(DataRate(min_rate)));
			rdmaHw->SetAttribute("Mtu", UintegerValue(packet_payload_size));
			rdmaHw->SetAttribute("MiThresh", UintegerValue(mi_thresh));
			rdmaHw->SetAttribute("VarWin", BooleanValue(var_win));
			rdmaHw->SetAttribute("FastReact", BooleanValue(fast_react));
			rdmaHw->SetAttribute("MultiRate", BooleanValue(multi_rate));
			rdmaHw->SetAttribute("SampleFeedback", BooleanValue(sample_feedback));
			rdmaHw->SetAttribute("TargetUtil", DoubleValue(u_target));
			rdmaHw->SetAttribute("RateBound", BooleanValue(rate_bound));
			rdmaHw->SetAttribute("DctcpRateAI", DataRateValue(DataRate(dctcp_rate_ai)));
			rdmaHw->SetAttribute("PowerTCPEnabled", BooleanValue(wien));
			rdmaHw->SetAttribute("PowerTCPdelay", BooleanValue(delayWien));
			rdmaHw->SetPintSmplThresh(pint_prob);
			// create and install RdmaDriver
			Ptr<RdmaDriver> rdma = CreateObject<RdmaDriver>();
			Ptr<Node> node = n.Get(i);
			rdma->SetNode(node);
			rdma->SetRdmaHw(rdmaHw);

			node->AggregateObject (rdma);
			rdma->Init();
			rdma->TraceConnectWithoutContext("QpComplete", MakeBoundCallback (qp_finish, fct_output));
		}
	}
	//cout << "ramdHW done.\n";

#endif
    // set ACK priority on hosts
	if (ack_high_prio)
		RdmaEgressQueue::ack_q_idx = 0;
	else
		RdmaEgressQueue::ack_q_idx = 3;

	// setup routing
	CalculateRoutes(n);
	cout << "Calculate Routes done.\n";
	SetRoutingEntries();
	cout << "Set Routing Entreis done.\n";

	//SetOpticalRoutes();
	//OpticalSwitch->ShowTableEntry();
    //
	// get BDP and delay 	获取server节点之间的BDP和延迟
	//
    maxRtt = maxBdp = 0;
	uint64_t minRtt = 1e9;  
    for (uint32_t i = 0; i < node_num; i++) {
		if (n.Get(i)->GetNodeType() != 0)	//switch
			continue;
		for (uint32_t j = 0; j < node_num; j++) {
			if (n.Get(j)->GetNodeType() != 0)
				continue;
			if (i == j)
				continue;
			uint64_t delay = pairDelay[n.Get(i)][n.Get(j)];
			uint64_t txDelay = pairTxDelay[n.Get(i)][n.Get(j)];
			uint64_t rtt = delay * 2 + txDelay;
			uint64_t bw = pairBw[i][j];
			uint64_t bdp = rtt * bw / 1000000000 / 8;
			pairBdp[n.Get(i)][n.Get(j)] = bdp;
			pairRtt[i][j] = rtt;
			if (bdp > maxBdp)
				maxBdp = bdp;
			if (rtt > maxRtt)
				maxRtt = rtt;
			if (rtt < minRtt)
				minRtt = rtt;
		}
	}
    printf("maxRtt=%lu maxBdp=%lu minRtt=%lu\n", maxRtt, maxBdp, uint64_t(minRtt));

    for (uint32_t i = 0; i < node_num; i++) {
		if (n.Get(i)->GetNodeType() && n.Get(i)->GetNodeType() != 3) { // switch
            if (n.Get(i)->GetNodeType() == 1){  //TorSwitch
                Ptr<TORSwitchNode> sw = DynamicCast<TORSwitchNode>(n.Get(i));
                sw->SetAttribute("CcMode", UintegerValue(cc_mode));
			    sw->SetAttribute("MaxRtt", UintegerValue(maxRtt));
			    sw->SetAttribute("PowerEnabled", BooleanValue(wien));
            }
            else{
                Ptr<SwitchNode> sw = DynamicCast<SwitchNode>(n.Get(i));
			    sw->SetAttribute("CcMode", UintegerValue(cc_mode));
			    sw->SetAttribute("MaxRtt", UintegerValue(maxRtt));
			    sw->SetAttribute("PowerEnabled", BooleanValue(wien));
            }
		}
	}

    Ipv4GlobalRoutingHelper::PopulateRoutingTables();

    NS_LOG_INFO("Create Applications.");

	Time interPacketInterval = Seconds(0.0000005 / 2);

    // maintain port number for each host
	for (uint32_t i = 0; i < node_num; i++) {
		if (n.Get(i)->GetNodeType() == 0)
			for (uint32_t j = 0; j < node_num; j++) {
				if (n.Get(j)->GetNodeType() == 0)
					portNumder[i][j] = 10000; // each host pair use port number from 10000
			}
	}

    flow_input.idx = 0;
	if (flow_num > 0) {
		ReadFlowInput();
		std::cout << flow_input.start_time << std::endl;
		Simulator::Schedule(Seconds(flow_input.start_time) - Simulator::Now(), ScheduleFlowInputs);
	}

    topof.close();
	tracef.close();
	double delay = 1.5 * minRtt * 1e-9; // 10 micro seconds
	Simulator::Schedule(Seconds(delay), PrintResults, switchDown, 6, delay);

    /*
    Here to set up the optical swtich
    Prepare to calculate demand matrix
    */
	// 设置 OCS 的属性
    OpticalSwitch->SetAttribute("ReconfigureDelay", UintegerValue(optical_reconfigure_delay));
    OpticalSwitch->SetAttribute("OswitchRate", UintegerValue(optical_oswitch_rate));
    OpticalSwitch->SetAttribute("PswitchRate", UintegerValue(optical_pswitch_rate));
	OpticalSwitch->SetAttribute("DurationExtentionRate", DoubleValue(optical_theta));
	if (MAN_OPTICAL == 1){ // 手动设置光切换矩阵
		OpticalSwitch->DemandMSetup(tors);
		std::ifstream opticalf;
		OpticalSwitch->SetDayorNight(false);
		opticalf.open(optical_conf_file);
		int tornum, dura_nums;
		double start_time;
		opticalf >> tornum >> dura_nums >> start_time;
		std::cout << "man: tornum " << tornum << " dura_nums " << dura_nums << " start_time " << start_time << '\n';
		uint64_t temp_dura_time;
		OpticalSwitch->man_init_configures(tornum, dura_nums);
		for (int i = 0; i < dura_nums; i++){
			opticalf >> temp_dura_time;
			OpticalSwitch->man_insert_duration_time(temp_dura_time);
		}
		bool temp_circuit;
		for (int i = 0; i < dura_nums; i++){
			for (int a = 0; a < tornum; a++){
				for (int b = 0; b < tornum; b++){
					opticalf >> temp_circuit;
					OpticalSwitch->man_set_configures(i, a, b, temp_circuit);
				}
			}
		}
		opticalf.close();
		OpticalSwitch->display_configures();
		Simulator::Schedule(Seconds(start_time), &OSwitchNode::Reconfigure, OpticalSwitch);
	}
	// 算法计算光切换矩阵 ***
	else{//auto calculate optical switch's demand and configures.
		/*Prepare the demand matrix of Optical switch.*/
		OpticalSwitch->DemandMSetup(tors);
		std::cout << "Going to set up the demand matrix.\n";
		FlowInput oflow_input = {0};
		std::ifstream oflowf;
		uint32_t oflow_num;
		oflowf.open(flow_file.c_str());
		oflowf >> oflow_num;
		double first_start;
		while (oflow_input.idx < oflow_num) {
			//std::cout << "input idx is " << oflow_input.idx <<" and oflow_num is " << oflow_num << "\n";
			oflowf >> oflow_input.src >> oflow_input.dst >> oflow_input.pg >> oflow_input.dport >> oflow_input.maxPacketCount >> oflow_input.start_time;
			std::cout << "Flow " << oflow_input.src << " " << oflow_input.dst << " " << oflow_input.pg << " " << oflow_input.dport << " " << oflow_input.maxPacketCount << " " << oflow_input.start_time << " " << Simulator::Now().GetSeconds() << std::endl;
			if(oflow_input.idx == 0){
				first_start = oflow_input.start_time;
			}
			Ptr<Node> srcToR, desToR;
			for (auto i = Tor2Server.begin(); i != Tor2Server.end(); i++){
				for (auto j = i->second.begin(); j != i->second.end(); j++){
					auto servernode = *j;
					if(servernode->GetId() == oflow_input.src){
						srcToR = i->first;
					}
					else if (servernode->GetId() == oflow_input.dst){
						desToR = i->first;
					}
				}
			}
			/* Optical need ToR to Interfaces map here.*/
			int in_indx,out_indx;
			for (auto i = nbr2if[OpticalSwitch].begin(); i != nbr2if[OpticalSwitch].end(); i++){
				if(i->first == srcToR){
					in_indx = i->second.idx;
				}
				else if( i->first == desToR ){
					out_indx = i->second.idx;
				}
			}
			//cout << "in: " << in_indx << " out: " << out_indx << '\n';
			//std::cout << "Calling SetFlowDemand(" << in_indx << ", " << out_indx << ", " << oflow_input.maxPacketCount * packet_payload_size << ")" << '\n';
			
			//OpticalSwitch->SetFlowDemand(in_indx, out_indx, oflow_input.maxPacketCount * packet_payload_size);
			OpticalSwitch->SetFlowDemand(in_indx, out_indx, (oflow_input.maxPacketCount * 8 / 1000000000 )); // byte to gbit
			oflow_input.idx++;
		}
		std::cout << "Demand matrix set.\n";
		OpticalSwitch->DisplayDemandM();
		//std::cout << "Calling CalculateConfigures.\n";
		OpticalSwitch->CalculateConfigures();
		//std::cout << "CalculateConfigures over.\n";
		OpticalSwitch->display_configures();
		Simulator::Schedule(Seconds(first_start), &OSwitchNode::Reconfigure, OpticalSwitch);
	}

	//----------------------------------------------------going to simulation-------------------------------------------
	//cout << "Printing next hops.\n";
	//show_nexthops();
	//cout << "Printing nbr2if.\n";
	//PrintNbr2If();
    std::cout << "Running Simulation.\n";
	NS_LOG_INFO("Run Simulation.");
	Simulator::Stop(Seconds(simulator_stop_time));
	Simulator::Run();
	Simulator::Destroy();
	NS_LOG_INFO("Done.");
	endt = clock();
	std::cout << "Simulation time cost: " << (double)(endt - begint) / CLOCKS_PER_SEC << "\n";
    return 0;
}
