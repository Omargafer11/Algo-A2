// ArcadiaEngine.cpp - STUDENT TEMPLATE
// TODO: Implement all the functions below according to the assignment requirements

#include "ArcadiaEngine.h"
#include <algorithm>
#include <queue>
#include <numeric>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>

using namespace std;

// =========================================================
// PART A: DATA STRUCTURES (Concrete Implementations)
// =========================================================

// --- 1. PlayerTable (Double Hashing) ---

class ConcretePlayerTable : public PlayerTable {
private:
    // TODO: Define your data structures here
    // Hint: You'll need a hash table with double hashing collision resolution
    struct PlayerEntry 
    {
        int id;
        string name;
        bool occupied; // To track if the slot is occupied
        bool deleted; // To track if the slot was deleted
    };
    static const int TABLE_SIZE = 101;
    PlayerEntry table[TABLE_SIZE]; // each element initialized to {0, "", false, false}


    // h1 --> Divison Method Hash Function
    int h1(int key) {
        return key % TABLE_SIZE;
    }


    // h2 --> Secondary Hash Function
    int h2(int key) {
        return 7 - (key % 7);
    }
public:
    ConcretePlayerTable() {
        // TODO: Initialize your hash table
        for (int i = 0; i < TABLE_SIZE; i++) {
            table[i] = {0, "", false, false};
        }
    }

    void insert(int playerID, string name) override {
        // TODO: Implement double hashing insert
        // Remember to handle collisions using h1(key) + i * h2(key)
        int index1 = h1(playerID);
        int index2 = h2(playerID);


        for(int i = 0; i < TABLE_SIZE; i++) 
        {
            int index = (index1 + i * index2) % TABLE_SIZE;

            if (!table[index].occupied || table[index].deleted) // Found an empty or deleted slot
            {
                table[index] = {playerID, name, true, false};
                return;
            }
            if (table[index].occupied && table[index].id == playerID) // Update existing player
            {
                table[index].name = name;
                return;
            }

        }
    }

    //searching using double hashing
    string search(int playerID) override 
    {
        // TODO: Implement double hashing search
        // Return "" if player not found
        int index1 = h1(playerID);
        int index2 = h2(playerID);

        for(int i = 0; i < TABLE_SIZE; i++) 
        {
            int index = (index1 + i * index2) % TABLE_SIZE;

            if (!table[index].occupied && !table[index].deleted) // Empty slot, player not found
            {
                return "";
            }
            if (table[index].occupied && table[index].id == playerID) // Found the player
            {
                return table[index].name;
            }
        }   
        return "";
    }
    void remove(int playerID) 
    {
        int index1 = h1(playerID);
        int index2 = h2(playerID);

        for(int i = 0; i < TABLE_SIZE; i++) 
        {
            int index = (index1 + i * index2) % TABLE_SIZE;

            if (!table[index].occupied && !table[index].deleted) // Empty slot, player not found
            {
                return;
            }
            if (table[index].occupied && table[index].id == playerID) // Found the player
            {
                table[index].occupied = false;// Mark as unoccupied
                table[index].deleted = true; // Mark as deleted
                return;
            }
        }   
    }
};

// --- 2. Leaderboard (Skip List) ---

class ConcreteLeaderboard  : public Leaderboard{
private:
    struct Node {
        int playerID;
        int score;
        vector<Node*> forward;
    };

    static const int MAX_LEVEL = 16;
    float P = 0.5;
    Node* header;
    int level;
    unordered_map<int, int> playerScores;

    int randomLevel() {
        int lvl = 1;
        while ((double)rand() / RAND_MAX < P && lvl < MAX_LEVEL)
            lvl++;
        return lvl;
    }

public:
    ConcreteLeaderboard() {
        level = 1;
        header = new Node{-1, INT_MAX, vector<Node*>(MAX_LEVEL, nullptr)};
    }

    void addScore(int playerID, int score) override{
        if (playerScores.count(playerID)) {
            removePlayer(playerID);
        }
        
        playerScores[playerID] = score;
        vector<Node*> update(MAX_LEVEL, nullptr);
        Node* cur = header;

        for (int i = level - 1; i >= 0; i--) {
            while (cur->forward[i] && (cur->forward[i]->score > score || 
                  (cur->forward[i]->score == score && cur->forward[i]->playerID < playerID))) {
                cur = cur->forward[i];
            }
            update[i] = cur;
        }

        int lvl = randomLevel();
        if (lvl > level) {
            for (int i = level; i < lvl; i++) update[i] = header;
            level = lvl;
        }

        Node* newNode = new Node{playerID, score, vector<Node*>(lvl, nullptr)};
        for (int i = 0; i < lvl; i++) {
            newNode->forward[i] = update[i]->forward[i];
            update[i]->forward[i] = newNode;
        }
    }

    void removePlayer(int playerID) override {
        if (playerScores.find(playerID) == playerScores.end()) return;

        int score = playerScores[playerID];
        vector<Node*> update(MAX_LEVEL, nullptr);
        Node* cur = header;

        for (int i = level - 1; i >= 0; i--) {
            while (cur->forward[i] && (cur->forward[i]->score > score || 
                  (cur->forward[i]->score == score && cur->forward[i]->playerID < playerID))) {
                cur = cur->forward[i];
            }
            update[i] = cur;
        }

        Node* target = cur->forward[0];
        if (target && target->playerID == playerID) {
            for (int i = 0; i < level; i++) {
                if (update[i]->forward[i] != target) break;
                update[i]->forward[i] = target->forward[i];
            }
            delete target;
            playerScores.erase(playerID);
        }

        while (level > 1 && header->forward[level - 1] == nullptr)
            level--;
    }

    vector<int> getTopN(int n)override {
        vector<int> result;
        Node* cur = header->forward[0];
        while (cur && n--) {
            result.push_back(cur->playerID);
            cur = cur->forward[0];
        }
        return result;
    }
};
// --- 3. AuctionTree (Red-Black Tree) ---

// Concrete implementation of AuctionTree using Red-Black Tree
class ConcreteAuctionTree : public AuctionTree {
private:
    // Node color definition
    enum Color { RED, BLACK };

    // Red-Black Tree node structure
    struct RBNode {
        int itemID;        // Key (used for BST ordering)
        int price;         // Associated item price
        Color color;       // Node color (RED or BLACK)
        RBNode *left;      // Pointer to left child
        RBNode *right;     // Pointer to right child
        RBNode *parent;    // Pointer to parent node

        // New nodes are always inserted as RED
        RBNode(int id, int p)
            : itemID(id), price(p), color(RED),
              left(nullptr), right(nullptr), parent(nullptr) {}
    };

    RBNode* root; // Root of the Red-Black Tree

    // Performs a left rotation around node x
    void leftRotate(RBNode* x) {
        RBNode* y = x->right;          // y becomes new parent of x
        x->right = y->left;

        if (y->left)
            y->left->parent = x;

        y->parent = x->parent;

        if (!x->parent)
            root = y;                 // x was root
        else if (x == x->parent->left)
            x->parent->left = y;
        else
            x->parent->right = y;

        y->left = x;
        x->parent = y;
    }

    // Performs a right rotation around node y
    void rightRotate(RBNode* y) {
        RBNode* x = y->left;           // x becomes new parent of y
        y->left = x->right;

        if (x->right)
            x->right->parent = y;

        x->parent = y->parent;

        if (!y->parent)
            root = x;                 // y was root
        else if (y == y->parent->left)
            y->parent->left = x;
        else
            y->parent->right = x;

        x->right = y;
        y->parent = x;
    }

    // Restores Red-Black Tree properties after insertion
    void fixInsert(RBNode* z) {
        // Continue while parent is RED (violation)
        while (z != root && z->parent->color == RED) {
            RBNode* p = z->parent;
            RBNode* g = p->parent;

            // Parent is left child of grandparent
            if (p == g->left) {
                RBNode* u = g->right;  // Uncle

                // Case 1: Uncle is RED (recolor)
                if (u && u->color == RED) {
                    p->color = u->color = BLACK;
                    g->color = RED;
                    z = g;
                } else {
                    // Case 2: Left-Right
                    if (z == p->right) {
                        z = p;
                        leftRotate(z);
                    }
                    // Case 3: Left-Left
                    p->color = BLACK;
                    g->color = RED;
                    rightRotate(g);
                }
            }
            // Parent is right child of grandparent (mirror cases)
            else {
                RBNode* u = g->left;   // Uncle

                if (u && u->color == RED) {
                    p->color = u->color = BLACK;
                    g->color = RED;
                    z = g;
                } else {
                    if (z == p->left) {
                        z = p;
                        rightRotate(z);
                    }
                    p->color = BLACK;
                    g->color = RED;
                    leftRotate(g);
                }
            }
        }
        // Root must always be BLACK
        root->color = BLACK;
    }

    // Returns the minimum node in a subtree (left-most node)
    RBNode* minimum(RBNode* n) {
        while (n->left)
            n = n->left;
        return n;
    }

    // Replaces subtree rooted at u with subtree rooted at v
    void transplant(RBNode* u, RBNode* v) {
        if (!u->parent)
            root = v;
        else if (u == u->parent->left)
            u->parent->left = v;
        else
            u->parent->right = v;

        if (v)
            v->parent = u->parent;
    }

    // Restores Red-Black Tree properties after deletion
    void fixDelete(RBNode* x, RBNode* parent) {
        while (x != root && (!x || x->color == BLACK)) {

            // x is left child
            if (x == parent->left) {
                RBNode* w = parent->right; // Sibling

                // Case 1: Sibling is RED
                if (w && w->color == RED) {
                    w->color = BLACK;
                    parent->color = RED;
                    leftRotate(parent);
                    w = parent->right;
                }

                // Case 2: Both sibling's children are BLACK
                if ((!w->left || w->left->color == BLACK) &&
                    (!w->right || w->right->color == BLACK)) {
                    if (w) w->color = RED;
                    x = parent;
                    parent = x->parent;
                }
                // Case 3 & 4
                else {
                    if (!w->right || w->right->color == BLACK) {
                        if (w->left) w->left->color = BLACK;
                        w->color = RED;
                        rightRotate(w);
                        w = parent->right;
                    }
                    w->color = parent->color;
                    parent->color = BLACK;
                    if (w->right) w->right->color = BLACK;
                    leftRotate(parent);
                    x = root;
                }
            }
            // x is right child (mirror cases)
            else {
                RBNode* w = parent->left;

                if (w && w->color == RED) {
                    w->color = BLACK;
                    parent->color = RED;
                    rightRotate(parent);
                    w = parent->left;
                }

                if ((!w->left || w->left->color == BLACK) &&
                    (!w->right || w->right->color == BLACK)) {
                    if (w) w->color = RED;
                    x = parent;
                    parent = x->parent;
                } else {
                    if (!w->left || w->left->color == BLACK) {
                        if (w->right) w->right->color = BLACK;
                        w->color = RED;
                        leftRotate(w);
                        w = parent->left;
                    }
                    w->color = parent->color;
                    parent->color = BLACK;
                    if (w->left) w->left->color = BLACK;
                    rightRotate(parent);
                    x = root;
                }
            }
        }
        if (x) x->color = BLACK;
    }

public:
    // Constructor initializes empty tree
    ConcreteAuctionTree() : root(nullptr) {}

    // Inserts or updates an item in the Red-Black Tree
    void insertItem(int itemID, int price) override {
        RBNode* z = new RBNode(itemID, price);
        RBNode* y = nullptr;
        RBNode* x = root;

        // Standard BST insertion
        while (x) {
            y = x;
            if (itemID < x->itemID) x = x->left;
            else if (itemID > x->itemID) x = x->right;
            else {
                // Item already exists → update price
                x->price = price;
                delete z;
                return;
            }
        }

        z->parent = y;
        if (!y) root = z;
        else if (itemID < y->itemID) y->left = z;
        else y->right = z;

        // Restore Red-Black Tree properties
        fixInsert(z);
    }

    // Deletes an item from the Red-Black Tree
    void deleteItem(int itemID) override {
        RBNode* z = root;

        // Search for the node
        while (z && z->itemID != itemID)
            z = (itemID < z->itemID) ? z->left : z->right;
        if (!z) return;

        RBNode* y = z;
        Color yColor = y->color;
        RBNode* x = nullptr;
        RBNode* xParent = nullptr;

        // Case 1: Node has at most one child
        if (!z->left) {
            x = z->right;
            xParent = z->parent;
            transplant(z, z->right);
        }
        else if (!z->right) {
            x = z->left;
            xParent = z->parent;
            transplant(z, z->left);
        }
        // Case 2: Node has two children
        else {
            y = minimum(z->right);     // In-order successor
            yColor = y->color;
            x = y->right;
            xParent = y->parent;

            transplant(y, y->right);
            y->right = z->right;
            y->right->parent = y;

            transplant(z, y);
            y->left = z->left;
            y->left->parent = y;
            y->color = z->color;
        }

        delete z;

        // Fix violations if a BLACK node was removed
        if (yColor == BLACK)
            fixDelete(x, xParent);
    }
};


// =========================================================
// PART B: INVENTORY SYSTEM (Dynamic Programming)
// =========================================================

int set_sum(const vector<int>& set) {
    int sum = 0;
    for (const int i : set)
        sum += i;
    return sum;
}

int InventorySystem::optimizeLootSplit(int n, vector<int>& coins) {
    // TODO: Implement partition problem using DP
    // Goal: Minimize |sum(subset1) - sum(subset2)|
    // Hint: Use subset sum DP to find closest sum to total/2
    int sum = set_sum(coins);
    vector<vector<int>> dp(coins.size() + 1, vector<int>(sum/2 + 1, 0));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= sum/2; j++) {
            if (coins[i - 1] > j) {
                dp[i][j] = dp[i - 1][j];
            }
            else {
                dp[i][j] = max(dp[i-1][j], coins[i-1] + dp[i-1][j - coins[i-1]]);
            }
        }
    }
    return abs(sum - 2*dp[n][sum/2]);
}

int InventorySystem::maximizeCarryValue(int capacity, vector<pair<int, int>>& items) {
    // TODO: Implement 0/1 Knapsack using DP
    // items = {weight, value} pairs
    // Return maximum value achievable within capacity
    int n = items.size();
    sort(items.begin(), items.end());
    vector<vector<int>> dp(n + 1, vector<int>(capacity + 1, 0));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= capacity; j++) {
            if (items[i - 1].first > j) {
                dp[i][j] = dp[i - 1][j];
            }
            else {
                dp[i][j] = max(dp[i-1][j], items[i-1].second + dp[i-1][j - items[i-1].first]);
            }
        }
    }
    return dp[n][capacity];
}

long long InventorySystem::countStringPossibilities(string s) {
    // TODO: Implement string decoding DP
    // Rules: "uu" can be decoded as "w" or "uu"
    //        "nn" can be decoded as "m" or "nn"
    // Count total possible decodings
    int n = s.length();
    long long MOD = 1e9 + 7;

    vector<long long> dp(n + 1);
    dp[0] = 1;
    dp[1] = 1;
    for (int i = 2; i <= n; i++) {
        dp[i] = (dp[i - 1] + dp[i - 2]) % MOD;
    }

    long long ans = 1;
    for (int i = 0; i < n; i++) {
        if (s[i] == 'u' || s[i] == 'n') {
            char current_char = s[i];
            int count = 0;
            while (i < n && s[i] == current_char) {
                count++;
                i++;
            }
            i--;

            ans = (ans * dp[count]) % MOD;
        }
    }

    return ans;
}

// =========================================================
// PART C: WORLD NAVIGATOR (Graphs)
// =========================================================


class Queue {
private:
    struct QueueNode {
        int data;
        QueueNode* next;
    };    
    QueueNode* front;
    QueueNode* rear;

public:
    Queue() {
        front = rear = nullptr;
    }

    bool isEmpty() {
        return front == nullptr;
    }

    void push(int x) {
        QueueNode* newNode = new QueueNode();
        newNode->data = x;
        newNode->next = nullptr;

        if (rear == nullptr) {
            front = rear = newNode;
            return;
        }

        rear->next = newNode;
        rear = newNode;
    }

    int pop() {
        if (isEmpty()) return -1;

        QueueNode* temp = front;
        int val = temp->data;

        front = front->next;
        if (front == nullptr)
            rear = nullptr;

        delete temp;
        return val;
    }
};

bool WorldNavigator::pathExists(int n, vector<vector<int>>& edges, int source, int dest) {

    // visited array
    bool visited[1000];
    for (int i = 0; i < n; i++)
        visited[i] = false;

    Queue q;

    q.push(source);
    visited[source] = true;

    while (!q.isEmpty()) {
        int current = q.pop();

        if (current == dest)
            return true;

        for (int i = 0; i < edges.size(); i++) {

            int neighbor = -1;

            if (edges[i][0] == current)
                neighbor = edges[i][1];
            else if (edges[i][1] == current)
                neighbor = edges[i][0];

            if (neighbor != -1 && !visited[neighbor]) {
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }

    return false;
}


int parent[1000];
int findParent(int x) {
    if (parent[x] == x)
        return x;
    return parent[x] = findParent(parent[x]);
}

void unionSet(int a, int b) {
    int pa = findParent(a);
    int pb = findParent(b);
    if (pa != pb)
        parent[pb] = pa;
}
struct Edge {
    int u, v;
    long long cost;
};

void sortEdges(Edge edges[], int m) {
    for (int i = 0; i < m - 1; i++) {
        for (int j = 0; j < m - i - 1; j++) {
            if (edges[j].cost > edges[j + 1].cost) {
                Edge temp = edges[j];
                edges[j] = edges[j + 1];
                edges[j + 1] = temp;
            }
        }
    }
}

long long WorldNavigator::minBribeCost(
    int n, int m, long long goldRate, long long silverRate,
    vector<vector<int>>& roadData) {

    Edge edges[1000];

    // حساب تكلفة كل طريق
    for (int i = 0; i < m; i++) {
        edges[i].u = roadData[i][0];
        edges[i].v = roadData[i][1];

        long long goldCost = roadData[i][2];
        long long silverCost = roadData[i][3];

        edges[i].cost = goldCost * goldRate + silverCost * silverRate;
    }

    sortEdges(edges, m);

    for (int i = 0; i < n; i++)
        parent[i] = i;

    long long totalCost = 0;
    int edgesUsed = 0;

    // Kruskal
    for (int i = 0; i < m; i++) {
        int u = edges[i].u;
        int v = edges[i].v;

        if (findParent(u) != findParent(v)) {
            unionSet(u, v);
            totalCost += edges[i].cost;
            edgesUsed++;
        }

        if (edgesUsed == n - 1)
            break;
    }

    if (edgesUsed != n - 1)
        return -1;

    return totalCost;
}

string WorldNavigator::sumMinDistancesBinary(int n, vector<vector<int>>& roads) {

    const long long INF = 1e18;
    long long dist[100][100];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                dist[i][j] = 0;
            else
                dist[i][j] = INF;
        }
    }

    for (int i = 0; i < roads.size(); i++) {
        int u = roads[i][0];
        int v = roads[i][1];
        int w = roads[i][2];

        if (w < dist[u][v]) {
            dist[u][v] = w;
            dist[v][u] = w;
        }
    }

    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }

    long long sum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (dist[i][j] != INF)
                sum += dist[i][j];
        }
    }

    if (sum == 0) return "0";

    string binary = "";
    while (sum > 0) {
        binary = char('0' + (sum % 2)) + binary;
        sum /= 2;
    }

    return binary;
}


// =========================================================
// PART D: SERVER KERNEL (Greedy)
// =========================================================

int ServerKernel::minIntervals(vector<char>& tasks, int n) {
    // TODO: Implement task scheduler with cooling time
    // Same task must wait 'n' intervals before running again
    // Return minimum total intervals needed (including idle time)
    // Hint: Use greedy approach with frequency counting
    int freq[26] = {0};
    for (char c : tasks) freq[c - 'A']++;

    // Find highest frequency
    int maxFreq = 0;
    for (int i = 0; i < 26; i++) if (freq[i] > maxFreq) maxFreq = freq[i];

    // Count how many tasks share this maxFreq
    int countMax = 0;
    for (int i = 0; i < 26; i++) if (freq[i] == maxFreq) countMax++;

    // Apply greedy formula
    int intervals = (maxFreq - 1) * (n + 1) + countMax;

    // The result must be at least the number of tasks
    return max(intervals, (int)tasks.size());
}

// =========================================================
// FACTORY FUNCTIONS (Required for Testing)
// =========================================================

extern "C" {
    PlayerTable* createPlayerTable() { 
        return new ConcretePlayerTable(); 
    }

    Leaderboard* createLeaderboard() { 
        return new ConcreteLeaderboard(); 
    }

    AuctionTree* createAuctionTree() { 
        return new ConcreteAuctionTree(); 
    }
}
