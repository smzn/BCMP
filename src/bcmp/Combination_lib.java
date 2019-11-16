package bcmp;

public class Combination_lib {
	private int n,c,index = 0,combi;
	private int[][] value;
	private int[] chosen;

	public Combination_lib(int n, int c) {
		this.n = n;
		this.c = c;
		//重複組合せの個数：cHn = c+n-1Cn = (c+n-1)!/(n!*(c-1)!)
		//combi = this.factorial(n + c - 1) / (this.factorial(n) * this.factorial(c - 1));
		//combi = this.getCombi_number2(c + n - 1, n);
		combi = this.getCombi_number3();
		value = new int[combi][c];
	}
	
	/*
	public int getCombi_number(int n, int r) {
		int sum1 = 0, sum2 = 0;
		for(int i = n - r + 1; i <= n; i++) {
			sum1 += Math.log(i);
		}
		for(int i = 2; i <= r; i++) {
			sum2 += Math.log(i);
		}
		double value = Math.exp(sum1 - sum2);
		System.out.println("対数重複組み合わせ" + value);
		return (int) Math.round(Math.exp(sum1 - sum2));
	}
	*/
	
	public int getCombi_number3() {
		int N = n + c - 1;
		int R = n;
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (int i = N - R + 1; i <= N; i++) {
			sum1 += Math.log(i);
		}
		for (int i = 2; i <= R; i++) {
			sum2 += Math.log(i);
		}
		double S = Math.exp(sum1 - sum2);
		// System.out.println(Math.round(S));
		return (int) Math.round(S);
	}
	
	public int getCombi_number2(int n, int r) {
		double sum = 1;
		for(int i = 1; i <= r; i++) {
			sum *= (n-r+i);
			sum /= i;
		}
		return (int) Math.round(sum);
	}
	
	public int factorial(int n){
		if(n == 0) return 1;
		return n * factorial(n-1);
	}
	
	public void GetRecursive(int n, int c_index, int pre_i) {
		if(c_index == c - 1) {
			value[index][c_index] = n;
			index ++;
			return;
		}
		else {
			for(int i = 0; i <= n; i++) {
				if(index == this.combi)break;
				if(c_index > 0) value[index][c_index - 1] = pre_i;
				value[index][c_index] = i;
				this.GetRecursive(n - i, c_index + 1, i);
				if(c_index == c-1 ) c_index = 0;
			}
		}
	}

	public int[][] getValue() {
		return value;
	}
	
	public int[][] getValueNc(int value[][], int nc[]){
		int nc_number = 0;
		//条件に当てはまる組合せを数える
		for(int i = 0; i < value.length; i++) {
			int flg = 1;
			for(int j = 0; j < c; j++) {
				if(value[i][j] > nc[j]) flg = -1;
			}
			if(flg == 1)nc_number++;
		}
		int valuenc[][] = new int[nc_number][c];
		int index = 0;
		for(int i = 0; i < value.length; i++) {
			int flg = 1;
			for(int j = 0; j < c; j++) {
				if(value[i][j] > nc[j]) flg = -1;
			}
			if(flg == 1) valuenc[index++] = value[i];
		}
		return valuenc;
	}
	
	
	/*ループでやる時
	public int[][] getCombination() {
		int value[][] = new int[6][c];
		int index = 0;
		for(int i = 0; i <= n; i++) {//第１層
			value[index][0] = i;
			for(int j = 0; j <= n-i;j++) {//第２層
				value[index][1] = j;
				//最下層
				value[index][2] = n - i - j;
				index++;
			}
		}
		return value;
	}
	*/
}
