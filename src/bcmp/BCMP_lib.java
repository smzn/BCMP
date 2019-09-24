package bcmp;

public class BCMP_lib {

	private int K,c,NK,N,n1,n2;
	private double lambda[], L[][], T[][], mu[][],rho[][];
	
	public BCMP_lib(int k, int c, int n, double mu[][]) {
		K = k; //ノード数
		this.c = c; //クラス数
		this.N = n; //系内人数
		NK = K * c -1;
		lambda = new double[c];
		L = new double[K][c];
		T = new double[K][c];
		this.mu = mu;
		rho = new double[c][K]; //利用率
		n1 = 60;//クラス別人数
		n2 = 40;
	}

	public double [] calcGauss(double[][] a, double[] b){
		int p;
		double pmax, s;
		double w[] = new double[NK];
		/* 前進消去（ピボット選択）*/
		for(int k = 0; k < NK-1; k++){  /* 第ｋステップ */
		      p = k;
		      pmax = Math.abs( a[k][k] );
		      for(int i = k+1; i < NK; i++){  /* ピボット選択 */
		         if(Math.abs( a[i][k] ) > pmax){
		            p = i;
		            pmax = Math.abs( a[i][k] );
		         }
		      }

		      if(p != k){  /* 第ｋ行と第ｐ行の交換　*/
		         for(int i = k; i < NK; i++){
		            /* 係数行列　*/
		            s = a[k][i];
		            a[k][i] = a[p][i];
		            a[p][i] = s;
		         }
		         /* 既知ベクトル */
		         s = b[k];
		         b[k] = b[p];
		         b[p] = s;
		      }
		/* 前進消去 */
		      for(int i = k +1; i < NK; i++){ /* 第ｉ行 */
		         w[i] = a[i][k] / a[k][k];
		         a[i][k] = 0.0;
		         /* 第ｋ行を-a[i][k]/a[k][k]倍して、第ｉ行に加える */
		         for(int j = k + 1; j < NK; j++){
		            a[i][j] = a[i][j] - a[k][j] * w[i];
		         }
		         b[i] = b[i] - b[k] * w[i];
		      }
		   }
		/* 後退代入 */
		      for(int i = NK - 1; i >= 0; i--){
		         for(int j = i + 1; j < NK; j++){
		            b[i] = b[i] - a[i][j] * b[j];
		            a[i][j] = 0.0;
		         }
		         b[i] = b[i] / a[i][i];
		         a[i][i] = 1.0;
		      }
		
		return b;
	}
	
	public void calcAverage(double alpha[]){
		int n = 0;
		while(n < N){
			n++;
			//Step2.1 Mean Response Time(滞在時間)
			for(int i = 0; i < K; i++){//Kはノード数
				for(int j = 0; j < c; j++) {//cはクラス数
					double sum = 0;
					for(int s = 0; s < c; s++) {//ノードiの全てのクラスのL[i][s]を求める
						sum += L[i][s]; 
					}
					T[i][j] = (sum + 1)/mu[i][j];
				}
			}
						
			//Step2.2 Throughput
			for(int i = 0; i < c;i++){
				double sum = 0;
				for(int j = 0; j < K; j++) {
					sum += alpha[i * K + j]*T[j][i];
				}
				if(i == 0) lambda[i] = n1/sum; //n1はクラス1の人数
				else lambda[i] = n2/sum; //n2はクラス2の人数
			}
						
			//Step2.3 系内人数
			for(int i = 0; i < K; i++){
				for(int j = 0; j < c; j++) {
					L[i][j] = lambda[j]*T[i][j]*alpha[i + j * K];
				}
			}
		}
	}

	public double[][] getT() {
		return T;
	}

	public double[] getLambda() {
		return lambda;
	}

	public double[][] getL() {
		return L;
	}

}
