
#include<iostream>
#include<vector>
#include "../alldef/matrix.cuh"
using namespace std;

void lower_upper_lemire(const vector<FLOAT >& a, int n, int r, vector<FLOAT >& l, vector<FLOAT >& u)
{
    int l_index = 0;
    int u_index = 0;
    int k = 2*r + 1;
    vector<int> q_u(n);
    vector<int> q_l(n);
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;
    for (int i = 0; i < k - 1 - r; i++) {
        while (head_l <= tail_l && a[q_l[tail_l]] >= a[i]) tail_l--;
        q_l[++tail_l] = i;
        while (head_u <= tail_u && a[q_u[tail_u]] <= a[i]) tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r ; i < n ; i++) {
        while (head_l <= tail_l && a[q_l[tail_l]] >= a[i]) tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k) head_l++;
        l[l_index++] = a[q_l[head_l]];

        while (head_u <= tail_u && a[q_u[tail_u]] <= a[i]) tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k) head_u++;
        u[u_index++] = a[q_u[head_u]];
    }

    for(int i = n;i < n + r;i++)
    {
        while (q_u[head_u] <= i - k) head_u++;
        u[u_index++] = a[q_u[head_u]];
        while (q_l[head_l] <= i - k ) head_l++;
        l[l_index++] = a[q_l[head_l]];
    }

}

void lower_upper_lemire(const FLOAT*  a, int n, int r, FLOAT*  l,FLOAT*  u)
{
    int l_index = 0;
    int u_index = 0;
    int k = 2*r + 1;
    vector<int> q_u(n);
    vector<int> q_l(n);
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;
    for (int i = 0; i < k - 1 - r; i++) {
        while (head_l <= tail_l && a[q_l[tail_l]] >= a[i]) tail_l--;
        q_l[++tail_l] = i;
        while (head_u <= tail_u && a[q_u[tail_u]] <= a[i]) tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r ; i < n ; i++) {
        while (head_l <= tail_l && a[q_l[tail_l]] >= a[i]) tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k) head_l++;
        l[l_index++] = a[q_l[head_l]];

        while (head_u <= tail_u && a[q_u[tail_u]] <= a[i]) tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k) head_u++;
        u[u_index++] = a[q_u[head_u]];
    }

    for(int i = n;i < n + r;i++)
    {
        while (q_u[head_u] <= i - k) head_u++;
        u[u_index++] = a[q_u[head_u]];
        while (q_l[head_l] <= i - k ) head_l++;
        l[l_index++] = a[q_l[head_l]];
    }

}