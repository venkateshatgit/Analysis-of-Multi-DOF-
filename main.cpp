#include<bits/stdc++.h>
using namespace std;

#define N 10
float stiff[N][N];
float mass[N][N];
float invOfMass[N][N];

void showmat(int n, float a[N][N], string x){
    cout<<"\n\t"<<x<<" = ";
    for(int i=0; i<n; ++i){
        cout<<"\n\t\t";
        for(int j=0; j<n; ++j){
            cout<<a[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<<"\n";
}

float determinant(float a[N][N], int n){

    float D=0;

    if(n==1)
        return a[0][0];

    else{
        float b[N][N];
        int s=1;

        for(int c=0; c<n; ++c){

            int p=0, q=0;
            for(int i=0; i<n; ++i){
                for(int j=0; j<n; ++j){
                    b[i][j]=0;

                    if(i!=0 && j!=c){
                        b[p][q]=a[i][j];
                        if(q<(n-2))
                            q++;
                        else{
                            q=0;
                            p++;
                        }
                    }
                }
            }
            D = D + s*a[0][c]*determinant(b, n-1);
            s=(-1)*s;
        }
    }
    return D;
}

int main(){


    while(true){



    int n, det=1;
    cout<<"\nEnter number of masses: ";
    cin>>n;
    float m[n];
    cout<<"\nEnter the value of masses:";
    for(int i=0; i<n; ++i){
        cin>>m[i];
        det = det*m[i];
    }



    for(int i=0; i<n; ++i){
        for(int j=0; j<n; ++j){
            if(i==j)
                mass[i][j]=m[i];
            else
                mass[i][j]=0;
        }
    }
    showmat(n, mass, "M");



    int k[n+1], p, d, s;

    cout<<"\nEnter number of supports: ";
    cin>>s;

    if(s==0){
        k[0]=0;
        k[n]=0;
        p=1;
        d=n-1;
    }

    else if(s==1){
        k[0]=0;
        p=1;
        d=n;
    }
    else{
        p=0;
        d=n+1;
    }


    cout<<"\nEnter "<<d<<" spring constants: ";

    if(p==0)
        d=n;

    for(int i=p; i<=d; ++i)
        cin>>k[i];

    for(int i=0; i<n+1; ++i)
        cout<<k[i]<<" ";

    for(int i=0; i<n; ++i){

        for(int j=0; j<n; ++j){

            if(i==j){
                stiff[i][j]=k[i]+k[i+1];

                if(i<n-1){
                    stiff[i+1][j]=(-1)*k[i+1];
                    stiff[i][j+1]=(-1)*k[i+1];
                }
            }
        }
    }

    cout<<"complete"<<endl;
    showmat(n, stiff, "K");


    for(int i=0; i<n; ++i){
        for(int j=0; j<n; ++j){
            if(i==j){
                invOfMass[i][j]=1/(m[i]);
            }
        }
    }

    showmat(n, invOfMass, "inv(M)");

    float mi[N][N];

    for(int i=0; i<n; ++i){

        for(int j=0; j<n; ++j){

            float sum=0;
            int a=i, b=0, c=0, d=j;
            for(int x=0; x<n; ++x){

                sum = sum + (invOfMass[a][b])*(stiff[c][d]);
                b++;
                c++;
            }

            mi[i][j]=sum;
        }
    }

    showmat(n, mi, "inv(M)*K");

    float solution[3];
    int t=0;

    for(float x=0; x<100; x+=0.1){
        float cop[N][N];

        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){
                cop[i][j]=mi[i][j];
            }
        }

        for(int i=0; i<n; ++i){
            for(int j=0; j<n; ++j){
                if(i==j){
                    cop[i][j]-=x;
                }
            }
        }

        float w=determinant(cop, n);

        if(w>=0 && w<1){
            solution[t]=x;
            t++;
            if(t==3)
                break;
        }

    }

    cout<<"\nEigen values are: ";

    if(t==0)
        cout<<"NO SOLUTION EXIST"<<endl;

    for(int i=0; i<t; ++i){

        cout<<solution[i]<<"\t";
        solution[i]=sqrt(solution[i]);
    }

    cout<<"\n\nCorresponding frequency is: ";

    for(int i=0; i<t; ++i){
        cout<<solution[i]<<"\t";
    }

    cout<<"\n\nRUN AGAIN (yes/no)"<<endl;

    string c;
    cin>>c;

    transform(c.begin(), c.end(), c.begin(), ::tolower);

    if(c=="no")
        break;

    }

return 0;
}
