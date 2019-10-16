package bcmp;

import java.util.Arrays;

public class MVA_lib {
	private int C,K,N,nclass[];//nclassはクラス毎の人数分布
	private double mu[][],lambda[][][],W[][][][],L[][][],alpha[][];//w,lの右側2つはクラス数で変動(課題)
	
	public MVA_lib(int c, int k, int n, int[] nclass, double[][] mu, double[][] alpha) {
		C = c; //クラス数
		K = k; //ノード数
		N = n; //客数
		this.nclass = nclass; //クラス毎の人数(最大数)
		this.mu = mu;
		this.alpha = alpha;
		lambda = new double[C][nclass[0]+1][nclass[1]+1];
		W = new double[c][k][nclass[0]+1][nclass[1]+1]; 
		L = new double[k][nclass[0]+1][nclass[1]+1];
	}
	
	public void getMVA() {
		for(int n = 1; n <= N; n++) {
			System.out.println("n = "+n);
			Combination_lib clib = new Combination_lib(n,C);
			clib.GetRecursive(n, 0, -1);
			int [][]value = clib.getValue();
			System.out.println("重複組合せ:" +Arrays.deepToString(value));
			int [][]valuenc = clib.getValueNc(value, nclass); //制限付きの総組合せ
			System.out.println("重複組合せ(制約付き):" +Arrays.deepToString(valuenc));
			for(int s = 0; s < valuenc.length; s++) { //nに対してのfeasible populationループ
				System.out.println("組合せ:" +Arrays.toString(valuenc[s]));
				//Step2 Wの計算
				for(int c = 0; c < C; c++) { //クラスのループ
					for(int k = 0; k < K; k++) { //ノードのループ
						if(c == 0) {
							if(valuenc[s][0] == 0) W[c][k][valuenc[s][0]][valuenc[s][1]] = 1/mu[c][k];
							else W[c][k][valuenc[s][0]][valuenc[s][1]] = 1/mu[c][k] * ( 1 + L[k][valuenc[s][0]-1][valuenc[s][1]]);
						}else if(c == 1){
							if(valuenc[s][1] == 0) W[c][k][valuenc[s][0]][valuenc[s][1]] = 1/mu[c][k];
							else W[c][k][valuenc[s][0]][valuenc[s][1]] = 1/mu[c][k] * ( 1 + L[k][valuenc[s][0]][valuenc[s][1]-1]);
						}
					}
				}// Step2終了
				//Step3 λの計算
				for(int c = 0; c < C; c++) {
					double sum = 0;
					for(int k = 0; k < K; k++) {
						sum += W[c][k][valuenc[s][0]][valuenc[s][1]] * alpha[c][k]; 
					}
					lambda[c][valuenc[s][0]][valuenc[s][1]] = nclass[c] / sum;
				}// Step3終了
				// Step4 Lの計算
				for(int k = 0; k < K; k++) {
					double sum = 0;
					for(int c = 0; c < C; c++) {
						sum += lambda[c][valuenc[s][0]][valuenc[s][1]] * W[c][k][valuenc[s][0]][valuenc[s][1]];
					}
					L[k][valuenc[s][0]][valuenc[s][1]] = sum; 
				}
			}
		}
	}

	public double[][][] getL() {
		return L;
	}

	public double[][][] getLambda() {
		return lambda;
	}

	public double[][][][] getW() {
		return W;
	}
	
}
