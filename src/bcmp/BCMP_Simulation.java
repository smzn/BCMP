package bcmp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class BCMP_Simulation {
	
	private double p[][], mu[][]; //muは1次元で渡す
	private int time;
	Random rnd = new Random();
	private int K; //ノード数
	private int N; //系内の客数
	private int C; //クラス数
	//ArrayList<Double> eventtime[]; 
	//ArrayList<String> event[];
	//ArrayList<Integer> queuelength[];
	//ArrayList<Integer> timequeue[];
	ArrayList<Integer> customer[];
	//private double timerate[][];
	//private double timerate2[][];
	//private double timeratecap[]; //キャパを超えている時間割合
	//private int capacity; //当面キャパは全ノードで一緒
	//private double maxlengthtime[];
	
	public BCMP_Simulation(double[][] p, double[][]mu, int time, int K, int N, int C) {
		this.p = p;
		this.mu = mu;
		this.time = time;
		this.K = K;
		this.N = N;
		this.C = C;
		//this.eventtime = new ArrayList[C*K];
		//this.event = new ArrayList[C*K];
		//this.queuelength = new ArrayList[C*K];
		//this.timequeue = new ArrayList[C*K];
		this.customer = new ArrayList[K];
		//for(int i = 0; i < eventtime.length; i++) eventtime[i] = new ArrayList<Double>();
		//for(int i = 0; i < event.length; i++) event[i] = new ArrayList<String>();
		//for(int i = 0; i < queuelength.length; i++) queuelength[i] = new ArrayList<Integer>();
		//for(int i = 0; i < timequeue.length; i++) timequeue[i] = new ArrayList<Integer>();
		for(int i = 0; i < customer.length; i++) customer[i] = new ArrayList<Integer>();
		
		//this.timerate = new double[C*K][N+1]; //0人の場合も入る
		//this.timerate2 = new double[N+1][C*K+1]; //0人からn人のn+1, 0拠点からk拠点でのk+1拠点
		//this.timeratecap = new double[C*K+1]; //timeratecap[i]キャパを超えたノード数がiの延べ時間、ノード数0-k(k+1個)
		//this.capacity = 2;
		//this.maxlengthtime = new double[C*K];
	}

	public double[][] getSimulation() {
		double service[] = new double[K];
		double result[][] = new double[2][K];
		int queue[] = new int[K]; //各ノードのサービス中を含むキューの長さ
		double elapse = 0;
		int error = 0, total = 0;
		
		//スタート時に2つのノードに分散
		for(int i = 0; i < this.N; i++) {
			//ノード00用
			//event[0].add("arrival");
			//queuelength[0].add(queue[0]);
			//eventtime[0].add(elapse); //(移動時間0)
			queue[0]++; //最初はノード0にn人いるとする
			//C=3の場合
			if(i%C == 0) customer[0].add(0); //クラス0として追加
			if(i%C == 1) customer[0].add(1); //クラス1として追加
			if(i%C == 2) customer[0].add(2);
			//event[1*K].add("arrival");
			//queuelength[1*K].add(queue[1*K]);
			//eventtime[1*K].add(elapse); //(移動時間0)
			//queue[1*K]++; //最初はノード0にn人いるとする
		}
		service[0] = this.getExponential(mu[0][0]); //先頭客のサービス時間設定
		//service[1*K] = this.getExponential(mu[1*K]); //先頭客のサービス時間設定
		double total_queue[] = new double[K]; //各ノードの延べ系内人数
		double total_queuelength[] = new double[K]; //待ち人数
		
		while(elapse < time) {
			double mini_service = 100000; //最小のサービス時間
			int mini_index = -1; //最小のサービス時間をもつノード
			int mini_class = -1; //サービス対象となる客のクラス
			for(int i = 0; i < K; i++) { //待ち人数がいる中で最小のサービス時間を持つノードを算出
				if( queue[i] > 0) {
					if( mini_service > service[i]) {
						mini_service = service[i];
						mini_index = i;
					}
				}
			}
			//System.out.println("mini_index = "+mini_index);
			//mini_indexの時のクラス
			mini_class = customer[mini_index].get(0);
			customer[mini_index].remove(0);//先頭を削除
			
			for(int i = 0; i < K; i++) { //ノードiから退去
				total_queue[i] += queue[i] * mini_service;
				if( queue[i] > 0) service[i] -= mini_service;
				if( queue[i] > 0 ) total_queuelength[i] += ( queue[i] - 1 ) * mini_service;
				else if ( queue[i] == 0 ) total_queuelength[i] += queue[i] * mini_service;
				//timerate[i][queue[i]] += mini_service;
			}
			
			//各ノードでの人数割合(同時滞在人数) 
			/*
			for(int i = 0; i < N+1; i++) { //0人からn人までのn+1
				int totalnumber = 0;
				for(int j = 0; j < queue.length; j++) { //0拠点からk拠点
					if(queue[j] == i) totalnumber ++;
				}
				timerate2[i][totalnumber] += mini_service;
			}
			//キャパを超えるノードの割合
			int totalnumber = 0; //ノード数
			for(int i = 0; i < C*K; i++) {
				if( queue[i] > capacity ) totalnumber ++;  
			}
			timeratecap[totalnumber] += mini_service;
			
			//イベント時間での待ち人数を登録
			//for(int i = 0; i < C*K; i++) timequeue[i].add(queue[i]);
			*/
			//event[mini_index].add("departure");
			//queuelength[mini_index].add(queue[mini_index]);
			queue[mini_index] --;
			elapse += mini_service;
			//eventtime[mini_index].add(elapse); //経過時間の登録はイベント後
			if( queue[mini_index] > 0) service[mini_index] = this.getExponential(mu[mini_class][mini_index]); 
			//退去後まだ待ち人数がある場合、サービス時間設定
			
			
			//ここから続きをやる(2019/11/16)
			//退去客の行き先決定
			double rand = rnd.nextDouble();
			//System.out.println("RND = "+rand);
			double sum_rand = 0;
			int destination_index = -1;
			for(int i = 0; i < p[0].length; i++) {
				sum_rand += p[mini_index][i];
				if( rand < sum_rand) {
					destination_index = i; //C*Kのどこか
					break;
				}
			}
			
			//行き先が決まらない場合
			if( destination_index == -1) {
				error ++;
				if(mini_index < K) {
					destination_index = K-1;//クラスの最後のノードにしておく
					//System.out.println("Class A Error");
				}
				else if (mini_index >= K ) {
					destination_index = K*2-1;
					//System.out.println("Class B Error");
				}
			}
			
			//destination_indexをクラスとノードに分ける
			int destination_class = destination_index / K; //クラス
			int destination_node = destination_index % K; //ノード
			
			// クラス内の移動を確認
			/*
			System.out.println("From : "+mini_index+" , To : "+destination_index);
			if(mini_index < K && destination_index < K)
				System.out.println("Class A");
			else if(mini_index >= K && destination_index >= K)
				System.out.println("Class B");
			else if(destination_index == -1) {
				System.out.println("Error");
			}*/
				
			//行先が決まらなかった場合(結構ある)
			//randが0.99以上になると行き先が決まらない
			/*
			total ++;
			if( destination_index == -1) {
				error ++;
				if(mini_index < K) {
					destination_index = K-1;//クラスの最後のノードにしておく
					//System.out.println("Class A Error");
				}
				else if (mini_index >= K ) {
					destination_index = K*2-1;
					//System.out.println("Class B Error");
				}
				//destination_index = p[0].length -1;
			}*/
			
			//event[destination_index].add("arrival");
			//queuelength[destination_index].add(queue[destination_index]);
			//eventtime[destination_index].add(elapse); //(移動時間0)
			//推移先で待っている客がいなければサービス時間設定(即時サービス)
			if(queue[destination_node] == 0) {
				service[destination_node] = this.getExponential(mu[destination_class][destination_node]);
			}
			queue[destination_node]++;	
			customer[destination_node].add(destination_class);
		} //Simulation roop終了
		
		for(int i = 0; i < K; i++) {
			result[0][i] = total_queue[i] / time; //平均系内人数
			result[1][i] = total_queuelength[i] / time; //平均待ち人数
		}
		//System.out.println("Error = "+error);
		return result;
		//System.out.println("Result:" +Arrays.deepToString(result));
	}
	
	//指数乱数発生
	public double getExponential(double param) {
		return - Math.log(1 - rnd.nextDouble()) / param;
	}
}
