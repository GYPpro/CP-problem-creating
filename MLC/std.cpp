
#include <iostream>
#include <queue>
#include <cstring>
using namespace std;

struct TrieNode {
    int cnt;
    TrieNode* children[26];
    TrieNode() : cnt(0) {
        memset(children, 0, sizeof(children));
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    TrieNode* root = new TrieNode();

    for (int i = 0; i < n; ++i) {
        string s;
        cin >> s;
        TrieNode* p = root;
        for (char c : s) {
            int idx = c - 'a';
            if (!p->children[idx]) {
                p->children[idx] = new TrieNode();
            }
            p = p->children[idx];
            p->cnt++;
        }
    }

    long long ans = 0;
    queue<TrieNode*> q;
    for (int i = 0; i < 26; ++i) {
        if (root->children[i]) {
            q.push(root->children[i]);
        }
    }

    while (!q.empty()) {
        TrieNode* node = q.front();
        q.pop();
        ans += 1LL * node->cnt * (node->cnt - 1) / 2;
        for (int i = 0; i < 26; ++i) {
            if (node->children[i]) {
                q.push(node->children[i]);
            }
        }
    }

    cout << ans << '\n';
    return 0;
}