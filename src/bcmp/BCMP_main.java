package bcmp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class BCMP_main {

	//課題1 クラス間移動がない場合、α11=1、α21=1と置かないといけない
	//今はalpha2として直接取り込んでいる
	//課題2 クラス数が2で固定
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int N = 10, K = 12, c = 2;
		int node_index[] = {14,15,17,18,21,24,27,29,31,34};
		int nc[] = {5,5};//各クラスの最大値
		double mu[][] = {{5,5,10,5,5,5,5,5,5,10,5,10}, {5,5,10,5,5,5,5,5,5,10,5,10}};//サービス率
		double mu_sim[] = {5,5,10,5,5,5,5,5,5,10,5,10,5,5,10,5,5,5,5,5,5,10,5,10};//シミュレーション用サービス率
		//double mu[][] = {{5,5},{5,5},{10,10},{5,5},{5,5},{5,5},{5,5},{5,5},{5,5},{10,10},{5,5},{10,10}};//サービス率
		double [][]r = new double[K * c][K * c];
		BCMP_main bmain = new BCMP_main();
		bmain.getCSV2("csv/transition.csv", K, c, r);
		System.out.println("推移確率行列" +Arrays.deepToString(r));
		
		double alpha[] = new double[K * c];//トラフィック方程式の解(α11=1とする)
		//トラフィック方程式を解く準備
		double alpha1[] = new double[K * c - 1];//トラフィック方程式α
		double r1[][] = new double[K * c -1][K * c -1];//Rを転置して対角要素を-1、１行と１列要素を削除
		double b1[] = new double[K * c -1]; //Rの１行目の2列目要素からK*c要素まで
		for(int i = 0; i < r.length -1; i++){
			for(int j = 0; j < r.length -1; j++){
				if( i == j ) {
					r1[i][j] = r[j + 1][i + 1] - 1;//転置して、対角要素は-1(１行、一列の各要素は除く) 
				}else {
					r1[i][j] = r[j + 1][i + 1];
				}
			}
		}
		for(int i = 0;i < r.length -1; i++){
			b1[i] = -r[0][i+1];
		}
		System.out.println("計算用行列左辺" +Arrays.deepToString(r1));
		System.out.println("計算用行列右辺" +Arrays.toString(b1));
		
		BCMP_lib blib = new BCMP_lib(K,c,N,mu);
		alpha1 = blib.calcGauss(r1, b1);//トラフィック方程式
		//α11=1としてalphaをalpha1から作る
		for(int i = 0 ; i < alpha.length; i++){
			if( i == 0) alpha[i] = 1;
			else alpha[i] = alpha1[i-1];
		}
		System.out.println("トラフィック方程式解α" +Arrays.toString(alpha));
		double alpha2[][] = {{1.0, 0.9760707373104454, 1.9469520689427025, 0.9502782204518215, 0.9765260877564881, 1.0196937083667734, 1.09062044987315, 1.0166426589311544, 0.9965026346573167, 1.702724752752357, 0.994438760077145, 1.4609981533055358},
				{1.0, 0.9760707373104454, 1.9469520689427025, 0.9502782204518215, 0.9765260877564881, 1.0196937083667734, 1.09062044987315, 1.0166426589311544, 0.9965026346573167, 1.702724752752357, 0.994438760077145, 1.4609981533055358}
		};
		
		Combination_lib clib = new Combination_lib(N,c);
		clib.GetRecursive(N,0, -1);
		int value[][] = clib.getValue();
		System.out.println("個数N = " +N);
		System.out.println("クラス数C = " +c);
		System.out.println("重複組合せ:個数" +value.length);
		System.out.println("重複組合せ:" +Arrays.deepToString(value));
		int valuenc[][] = clib.getValueNc(value, nc);
		System.out.println("制約(最大数):" +Arrays.toString(nc));
		System.out.println("重複組合せ(制約付き):個数" +valuenc.length);
		System.out.println("重複組合せ(制約付き):" +Arrays.deepToString(valuenc));
		
		MVA_lib mlib = new MVA_lib(c, K, N, nc, mu,alpha2);
		mlib.getMVA();
		double L[][][] = mlib.getL();
		double W[][][][] = mlib.getW();
		double lambda[][][] = mlib.getLambda();
		System.out.println("W:" +Arrays.deepToString(W));
		System.out.println("λ:" +Arrays.deepToString(lambda));
		System.out.println("L:" +Arrays.deepToString(L));
		
		for(int i = 0; i < W.length; i++) {
			for(int j = 0; j < W[i].length; j++) {
				System.out.println("W["+i+"]["+j+"]["+nc[0]+"]["+nc[1]+"]= "+W[i][j][nc[0]][nc[1]]);
			}
		}
		for(int i = 0; i < lambda.length; i++) {
			System.out.println("lambda["+i+"]["+nc[0]+"]["+nc[1]+"]= "+lambda[i][nc[0]][nc[1]]);
		}
		for(int i = 0; i < L.length; i++) {
			System.out.println("L["+i+"]["+nc[0]+"]["+nc[1]+"]= "+L[i][nc[0]][nc[1]]);
		}
		
		//Simulation
		BCMP_Simulation slib = new BCMP_Simulation(r, mu_sim, 300000, node_index, K, N, c);
		double [][] result = slib.getSimulation();
		System.out.println("Result:" +Arrays.deepToString(result));
	}

	public void getCSV2(String path, int K, int c, double r[][]) {
		//CSVから取り込み
		try {
			File f = new File(path);
			BufferedReader br = new BufferedReader(new FileReader(f));
				 
			String[][] data = new String[K * c][K * c]; 
			String line = br.readLine();
			for (int i = 0; line != null; i++) {
				data[i] = line.split(",", 0);
				line = br.readLine();
			}
			br.close();
			// CSVから読み込んだ配列の中身を処理
			for(int i = 0; i < data.length; i++) {
				for(int j = 0; j < data[0].length; j++) {
					r[i][j] = Double.parseDouble(data[i][j]);
				}
			} 
		} catch (IOException e) {
			System.out.println(e);
		}
		//CSVから取り込みここまで	
	}
}
