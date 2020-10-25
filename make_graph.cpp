#include <bits/stdc++.h>
using namespace std;

int main()
{
    int x, y;
    int n=5;
    int m=5;
    int cur = 0;

    freopen("input.txt", "w", stdout);
    cout << n << " " << m << endl;

    vector<int> *v = new vector<int>[n+1];
    for (int i = 0; i < m; i++) {
        do {
            x = rand() % n;
            y = rand() % n;
        } while (x == y);

        v[x].push_back(y);
        v[y].push_back(x);
    }

    for(int i=0; i<n; i++)
        for(int node : v[i])
            cout<<node<<" ";
    cout<<"\n";

    cout<<"0 ";
    for(int i=0; i<n; i++) {
        cur += v[i].size();
        cout<<cur<<" ";
    }
    cout<<"\n";
}