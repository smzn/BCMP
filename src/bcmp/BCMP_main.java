package bcmp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class BCMP_main {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int N = 100, K = 12, c = 2;
		//double mu[][] = {{5,5,10,5,5,5,5,5,5,10,5,10}, {5,5,10,5,5,5,5,5,5,10,5,10}};//サービス率
		double mu[][] = {{5,5},{5,5},{10,10},{5,5},{5,5},{5,5},{5,5},{5,5},{5,5},{10,10},{5,5},{10,10}};//サービス率
		double [][]r = new double[K * c][K * c];
		BCMP_main bmain = new BCMP_main();
		bmain.getCSV2("csv/transition_class.csv", K, c, r);
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
		
		blib.calcAverage(alpha);
		double T[][] = blib.getT();
		double lambda[] = blib.getLambda();
		double L[][] = blib.getL();
		System.out.println("平均滞在時間" +Arrays.deepToString(T));
		System.out.println("Throughput" +Arrays.toString(lambda));
		System.out.println("平均系内人数" +Arrays.deepToString(L));
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
