#include <iostream>
#include <vector>
#include <cstring>

using namespace std;

bool bpm(int u, vector<vector<int>>& bpGraph, vector<bool>& seen, vector<int>& matchR) {
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

int maxBPM(vector<vector<int>>& bpGraph, vector<int>& matchR) {
    int U = bpGraph.size();
    int V = bpGraph[0].size();
    for (int u = 0; u < U; u++) {
        vector<bool> seen(V, false);
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

bool isPerfectMatch(vector<int>& matchR, int U, int V) {
    int matchedCount = 0;
    for (int v = 0; v < V; v++) {
        if (matchR[v] != -1) {
            matchedCount++;
        }
    }
    return matchedCount == U && matchedCount == V;
}

void printResultMatrix(vector<int>& matchR, int U, int V) {
    vector<vector<int>> result(U, vector<int>(V, 0));
    for (int v = 0; v < V; v++) {
        if (matchR[v] != -1) {
            result[matchR[v]][v] = 1;
        }
    }
    
    cout << "Result Matrix:" << endl;
    for (int i = 0; i < U; i++) {
        for (int j = 0; j < V; j++) {
            cout << result[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    int U = 5; // number of vertices on the left side
    int V = 5; // number of vertices on the right side

    int matrix[5][5] = { {0, 0, 0, 1, 0},
                         {0, 0, 0, 0, 1},
                         {1, 0, 0, 0, 0},
                         {0, 1, 0, 0, 0}, 
                         {0, 0, 1, 0, 0} };

    vector<vector<int>> bpGraph(U, vector<int>(V, 0));

    for (int i = 0; i < U; i++) {
        for (int j = 0; j < V; j++) {
            bpGraph[i][j] = matrix[i][j];
        }
    }

    vector<int> matchR(V, -1);

    int maxMatch = maxBPM(bpGraph, matchR);

    cout << "Maximum number of matches: " << maxMatch << endl;

    printResultMatrix(matchR, U, V);

    if (isPerfectMatch(matchR, U, V)) {
        cout << "The maximum matching is a perfect matching." << endl;
    } else {
        cout << "The maximum matching is not a perfect matching." << endl;
    }

    return 0;
}
