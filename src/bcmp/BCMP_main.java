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
	//20191118
	//シミュレーションを変更してだいぶ値が近くなったがもう少し何とかしたい
	//クラス間移動を許した時、平均系内人数がL[6]=29->L[6]=27に減少する。
	//クラス間移動なしと移動あり(クラスだけ変えて、移動予定ノードは変わらない)で各クラスの内容が同じなので、同じ値になるはずだけど、ならない
	//シミュレーションは同じ値になる
	//トラフィック方程式の部分が怪しいと思うのでチェックする？
	//20191117_計算結果.xlsxに記載
	//MVAは集中する拠点があると値の変化が大きい(調査した方がいい)
	//クラス間移動なしでやった場合、拠点での系内人数は変化しないはずだけど、MVAでは変化してしまった(集中する場合)
	//閉鎖型の場合、α11=1と置いていて、これは客がノード1を訪問した場合に拠点iを何度訪れたかの平均訪問回数αciなので、クラス間移動がない場合は概念的にはおかしいけど、α21=1とおいて、クラス2でのαを作ってやる
	//ならした値でやればMVAでも良さそう。
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int N = 100, K = 24, c = 2;
		int nc[] = {50,50};//各クラスの最大値
		//12拠点
		//double mu[][] = {{5,5,10,5,5,5,7,5,5,10,5,10}, {5,5,10,5,5,5,7,5,5,10,5,10}};//サービス率
		//24拠点
		double mu[][] = {{5,5,10,5,5,5,7,5,5,10,5,10,5,5,10,5,5,5,7,5,5,10,5,10}, {5,5,10,5,5,5,7,5,5,10,5,10,5,5,10,5,5,5,7,5,5,10,5,10}};
		double [][]r = new double[K * c][K * c];
		double alpha[] = new double[K * c];//トラフィック方程式の解(α11=1とする)
		double alpha2[][] = new double[c][K];//クラス別2次元配列
		
		//(1) 推移確率行列の取り込み
		BCMP_main bmain = new BCMP_main();
		bmain.getCSV2("csv/transition2_24.csv", K, c, r);
		System.out.println("推移確率行列" +Arrays.deepToString(r));
		
		//(2)トラフィック方程式を解く
		double alpha1[] = new double[K * c - 1];//トラフィック方程式α
		double r1[][] = new double[K * c -1][K * c -1];//Rを転置して対角要素を-1、１行と１列要素を削除
		double b1[] = new double[K * c -1]; //Rの１行目の2列目要素からK*c要素まで
		BCMP_lib blib = new BCMP_lib(K,c,N,mu);
		blib.getPretrafiic(r);//行列の変形
		r1 = blib.getR1();
		b1 = blib.getB1();
		alpha1 = blib.calcGauss(r1, b1);//トラフィック方程式
		blib.getAlpha_value(alpha1);//alpha1を元の形に
		alpha = blib.getAlpha();//alpha1次元
		alpha2 = blib.getAlpha_class();//alphaクラス別
		System.out.println("トラフィック方程式解α(クラス別)" +Arrays.deepToString(alpha2));
		
		//α11=1としてalphaをalpha1から作る
		for(int i = 0 ; i < alpha.length; i++){
			if( i == 0) alpha[i] = 1;
			else alpha[i] = alpha1[i-1];
		}
		System.out.println("トラフィック方程式解α" +Arrays.toString(alpha));
		//double alpha2[][] = new double[c][K];
		for(int i = 0; i < c; i++) {
			for(int j = 0; j < K; j++) {
				alpha2[i][j] = alpha[i * K + j];
			}
		}
		//System.out.println("トラフィック方程式解α2" +Arrays.deepToString(alpha2));
		/*
		double alpha3[][] = {{1.0, 0.9760707373104454, 1.9469520689427025, 0.9502782204518215, 0.9765260877564881, 1.0196937083667734, 1.09062044987315, 1.0166426589311544, 0.9965026346573167, 1.702724752752357, 0.994438760077145, 1.4609981533055358},
				{1.0, 0.9760707373104454, 1.9469520689427025, 0.9502782204518215, 0.9765260877564881, 1.0196937083667734, 1.09062044987315, 1.0166426589311544, 0.9965026346573167, 1.702724752752357, 0.994438760077145, 1.4609981533055358}
		};
		*/
		
		//組み合わせ計算の階乗をやらない計算をする(Nを大きくできるように)
		//http://slpr.sakura.ne.jp/qp/large-number/
		/*
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
		*/
		
		//(3) MVAでの計算
		MVA_lib mlib = new MVA_lib(c, K, N, nc, mu,alpha2);
		mlib.getMVA();
		double L[][][] = mlib.getL();
		double W[][][][] = mlib.getW();
		double lambda[][][] = mlib.getLambda();
		double L_node[] = new double[K];//平均系内人数最終結果
		System.out.println("W:" +Arrays.deepToString(W));
		System.out.println("λ:" +Arrays.deepToString(lambda));
		System.out.println("L[6]:" +Arrays.deepToString(L[6]));
		
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
			L_node[i] = L[i][nc[0]][nc[1]];
		}
		
		//(4) Simulation
		int time = 300000;
		BCMP_Simulation slib = new BCMP_Simulation(r, mu, time, K, N, c);
		System.out.println("Simulation Start");
		double [][] result = slib.getSimulation();
		System.out.println("Result:" +Arrays.deepToString(result));
		
		//(5) DBへ格納
		int combination_id = 2, transition_id = 1;
		MySQL mysql = new MySQL(combination_id, transition_id);
		//mysql.insertL(L_node);
		//mysql.insertSimulationL(result[0], time);
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
