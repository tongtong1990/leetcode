import java.awt.datatransfer.StringSelection;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Stack;

class ListNode{
	int val;
	ListNode next;
	ListNode(int x) {
		val = x;
		next = null;
	}
}

class TreeNode {
	int val;
	TreeNode left;
	TreeNode right;
	TreeNode(int x) { val = x; }
}

class Interval{
	int start;
	int end;
	Interval() { start = 0; end = 0; }
	Interval(int s, int e) { start = s; end = e; }
}

class TreeLinkNode{
	int val;
	TreeLinkNode left, right, next;
	TreeLinkNode(int x) { val = x; }
}

public class practice {
	//3sum
    public ArrayList<ArrayList<Integer>> threeSum(int[] num) {
        // Start typing your Java solution below
        // DO NOT write main() function
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        if(num.length < 3)	return result;
        Arrays.sort(num);
        int rest, prev=num[0], i=0;
        while(i<num.length-2){
        	if(num[i] > 0)
        		return result;
        	else{
        		rest = -num[i];
        		int left = i+1;
        		int right = num.length-1;
        		while(left < right){
        			int leftValue = num[left];
        			int rightValue = num[right];
        			if(leftValue + rightValue == rest){
        				ArrayList<Integer> al = new ArrayList<Integer>();
        				al.add(-rest);
        				al.add(num[left]);
        				al.add(num[right]);
        				result.add(al);
        				while(num[left] == leftValue && left < right)	left ++;
        				while(num[right] == rightValue && left < right) right --;
        			}
        			else if(num[left] + num[right] > rest)
        				right --;
        			else
        				left ++;
        		}
        	}
        	prev = num[i];
        	while(num[i] == prev && i<num.length-2)	i++;
        }
        return result;
    }
    
    //3sum closest
    public int threeSumClosest(int[] num, int target) {
        // Start typing your Java solution below
        // DO NOT write main() function
        int result = 0;
        if(num.length < 3)	return result;
        Arrays.sort(num);
        int d = Integer.MAX_VALUE;
        for(int i=0;i<num.length-2;i++){
        	int left = i+1;
        	int right = num.length-1;
        	int rest = target - num[i];
        	while(left < right){
        		if(num[left] + num[right] == rest)	return target;
        		if(Math.abs(rest-num[left]-num[right]) < d){
        			d = Math.abs(rest-num[left]-num[right]);
        			result = num[i] + num[left] + num[right];
        		}
        		if(num[left] + num[right] > rest)
        			right --;
        		else
        			left ++;
        	}
        }
        return result;
    }
    
    //4sum
    public ArrayList<ArrayList<Integer>> fourSum(int[] num, int target) {
        // Start typing your Java solution below
        // DO NOT write main() function
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        if(num.length < 4)	return result;
        Arrays.sort(num);
        int first = 0,second;
        while(first < num.length-3){
        	second = first + 1;
        	while(second < num.length-2){
            	int left = second + 1, right = num.length-1;
            	while(left < right){
            		int tmp_sum = num[left] + num[right] + num[first] + num[second];
            		if(tmp_sum == target){
            			ArrayList<Integer> al = new ArrayList<Integer>();
            			al.add(num[first]); al.add(num[second]); al.add(num[left]); al.add(num[right]);
            			result.add(al);
            			int leftV = num[left], rightV = num[right];
            			while(num[left] == leftV && left < right)
            				left ++;
            			while(num[right] == rightV && left < right)
            				right --;
            		}
            		else if(tmp_sum < target){
            			left ++;
            		}
            		else{
            			right --;
            		}
            	}
            	int secV = num[second];
            	while(num[second] == secV && second < num.length-2)
            		second ++;
        	}
        	int firstV = num[first];
        	while(num[first] == firstV && first < num.length-3)
        		first ++;
        }
        return result;
    }
    
    //add binary
    public String addBinary(String a, String b) {
        // Start typing your Java solution below
        // DO NOT write main() function
        String result = "";
        int carry = 0;
        char [] aa = a.toCharArray();
        char [] bb = b.toCharArray();
        if(aa.length < bb.length){
        	char []tmp = bb;
        	bb = aa;
        	aa = tmp;
        }
        int [] inta = new int[aa.length];
        int [] intb = new int[bb.length];
        for(int i=0;i<aa.length;i++)
        	inta[i] = aa[i] - '0';
        for(int i=0;i<bb.length;i++)
        	intb[i] = bb[i] - '0';
        int [] res = new int[aa.length];
        int i = 0;
        for(i=0;i<bb.length;i++){
        	int a_bit = inta[aa.length-1-i];
        	int b_bit = intb[bb.length-1-i];
        	if(carry == 0){
        		carry = a_bit & b_bit;
        		res[aa.length-1-i] = a_bit ^ b_bit;
        	}
        	else{
        		carry = a_bit | b_bit;
        		res[aa.length-1-i] = 1 - (a_bit ^ b_bit);
        	}
        }
        for(;i<aa.length;i++){
        	if(carry == 0){
        		carry = 0;
        		res[aa.length-1-i] = inta[aa.length-1-i];
        	}
        	else{
        		res[aa.length-1-i] = inta[aa.length-i-1] ^ carry;
        		carry = inta[aa.length-1-i] & carry;
        	}
        }
        char []res_char = new char[aa.length];
        for(i=0;i<aa.length;i++)
        	res_char[i] = (char)(res[i] + '0');
        result = String.valueOf(res_char);
        if(carry == 1)	result = '1' + result;
        return result;
    }
	
    //add two numbers
    public ListNode addTwoNumbers(ListNode l1, ListNode l2) {
        // Start typing your Java solution below
        // DO NOT write main() function
        ListNode result = null;
        ListNode prev = null;
        int carry = 0;
        while(l1 != null && l2 != null ){
        	ListNode tmp = new ListNode((l1.val + l2.val + carry) % 10);
        	carry = (l1.val + l2.val + carry) / 10;
        	if(prev == null){
        		result = prev = tmp;
        	}
        	else{
        		prev.next = tmp;
        		prev = tmp;
        	}
        	l1 = l1.next;
        	l2 = l2.next;
        }
        if(l2 != null)
        	l1 = l2;
        while(l1 != null){
        	ListNode tmp = new ListNode((l1.val + carry) % 10);
        	carry = (l1.val + carry) / 10;
        	if(prev == null)
        		result = prev = tmp;
        	else{
        		prev.next = tmp;
        		prev = tmp;
        	}
        	l1 = l1.next;
        }
        if(carry != 0){
        	ListNode tmp = new ListNode(carry);
        	prev.next = tmp;
        }
        return result;
    }
    
    //anagrams
    public ArrayList<String> anagrams(String[] strs) {
        // Start typing your Java solution below
        // DO NOT write main() function
        ArrayList<String> result = new ArrayList<String> ();
        HashSet<String> hs1 = new HashSet<String>();
        HashSet<String> hs2 = new HashSet<String>();
        String []strs2 = new String[strs.length];
        for(int i=0; i<strs.length; i++){
        	char []tmp = strs[i].toCharArray();
        	Arrays.sort(tmp);
        	String tmpS = new String(tmp);
        	strs2[i] = tmpS;
        	if(!hs1.contains(tmpS))
        		hs1.add(tmpS);
        	else if(!hs2.contains(tmpS))
        		hs2.add(tmpS);
        }
        for(int i=0; i<strs.length; i++){
        	if(hs2.contains(strs2[i]))
        		result.add(strs[i]);
        }
        return result;
    }
    
    //balanced binary tree
    public boolean isBalanced(TreeNode root) {
        // Start typing your Java solution below
        // DO NOT write main() function
        boolean result = false;
        if(isBalancedRec(root) == -1)
        	return false;
        return true;
    }
    public int isBalancedRec(TreeNode root){
    	if(root == null)	return 0;
    	int left = isBalancedRec(root.left);
    	int right = isBalancedRec(root.right);
    	if(left == -1 || right == -1 || left-right > 1 || right-left > 1)	return -1;
    	return left > right ? left+1 : right+1;
    }

    //best time to buy and sell stock
    public int maxProfit(int[] prices) {
        // Start typing your Java solution below
        // DO NOT write main() function
        int result = 0;
        if(prices.length < 2)	return result;
        int lowest = prices[0];
        for(int i=0; i<prices.length; i++){
        	if(prices[i] < lowest){
        		lowest = prices[i];
        	}
        	else if(prices[i] - lowest > result){
        		result = prices[i] - lowest;
        	}
        }
        return result;
    }

    //best time to buy and sell stock 2
    public int maxProfit2(int[] prices) {
        // Start typing your Java solution below
        // DO NOT write main() function
        int result = 0;
        if(prices.length < 2)	return result;
        int prev = prices[0];
        int lowest = prices[0];
        for(int i=0;i<prices.length;i++){
        	if(prices[i] > prev){
        		prev = prices[i];
        	}else{
        		result += (prev - lowest);
        		prev = lowest = prices[i];
        	}
        }
        result += (prev - lowest);
        return result;
    }
	
    //best time to buy and sell stock 3
    public int maxProfit3(int[] prices) {
    	//D&C & DP
    	//key : must be two transactions ( must not be single transaction. can prove it)
        int result = 0;
        if(prices.length < 2)	return result;
        int []left = new int[prices.length];
        int []right = new int[prices.length];
        int lowest = prices[0];
        int tmp_result = 0;
        for(int i=0;i<prices.length;i++){
        	tmp_result = (prices[i] - lowest) > tmp_result ? prices[i] - lowest : tmp_result;
        	left[i] = tmp_result;
        	if(prices[i] < lowest)
        		lowest = prices[i];
        }
        tmp_result = 0;
        int highest = prices[prices.length-1];
        for(int i=prices.length-1 ; i>=0 ; i--){
        	tmp_result = (highest - prices[i]) > tmp_result ? highest - prices[i] : tmp_result;
        	right[i] = tmp_result;
        	if(prices[i] > highest)
        		highest = prices[i];
        }
        for(int i=0;i<prices.length-1;i++){
        	if(left[i] + right[i] > result)
        		result = left[i] + right[i];	// instead of left[i] + right[i+1] because can sell and then buy, this is because why transcations must be two.
        }
        return result;
    }

    //binary tree in order traversal
    public ArrayList<Integer> inorderTraversal(TreeNode root) {
    	ArrayList<Integer> result = new ArrayList<Integer>();
    	if(root == null)	return result;
    	class tmp_node{
    		TreeNode node;
    		boolean isVisited = false;
    		tmp_node(TreeNode n) { node = n; }
    	}
    	tmp_node tn = new tmp_node(root);
    	Stack<tmp_node> st = new Stack<tmp_node>();
    	st.push(tn);
    	while(!st.isEmpty()){
    		tmp_node cur = st.peek();
    		if(cur.isVisited){
    			st.pop();
    			result.add(cur.node.val);
    			if(cur.node.right != null){
    				tmp_node new_node = new tmp_node(cur.node.right);
    				st.push(new_node);
    			}
    		}
    		else{
    			cur.isVisited = true;
    			if(cur.node.left != null){
    				tmp_node new_node = new tmp_node(cur.node.left);
    				st.push(new_node);
    			}
    		}
    	}
    	
    	return result;
    }
    
    //binary tree level order traversal 
    public ArrayList<ArrayList<Integer>> levelOrderBottom(TreeNode root) {
        // Start typing your Java solution below
        // DO NOT write main() function
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        Stack<ArrayList<Integer>> st = new Stack<ArrayList<Integer>>();
        Queue<TreeNode> q = new LinkedList<TreeNode>();
        int thisLayerNum = 1, nextLayerNum = 0;
        if(root == null)	return result;
        ArrayList<Integer> thisLayerArray = new ArrayList<Integer>();
        q.offer(root);
        while(!q.isEmpty()){
        	TreeNode cur = q.poll();
        	thisLayerArray.add(cur.val);
        	thisLayerNum --;
        	if(cur.left != null){
        		nextLayerNum ++;
        		q.offer(cur.left);
        	}
        	if(cur.right != null){
        		nextLayerNum ++;
        		q.offer(cur.right);
        	}
        	if(thisLayerNum == 0){
        		thisLayerNum = nextLayerNum;
        		nextLayerNum = 0;
        		st.push(thisLayerArray);
        		thisLayerArray = new ArrayList<Integer>();
        	}
        }
        while(!st.isEmpty())
        	result.add(st.pop());
        return result;
    }

    //binary tree level order traversal 2
    public ArrayList<ArrayList<Integer>> levelOrder(TreeNode root) {
    	ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
    	if(root == null)	return result;
    	int thisLayerNum = 1, nextLayerNum = 0;
    	Queue<TreeNode> q = new LinkedList<TreeNode>();
    	q.offer(root);
    	ArrayList<Integer> thisLayerInt = new ArrayList<Integer>();
    	while(!q.isEmpty()){
    		TreeNode cur = q.poll();
    		thisLayerInt.add(cur.val);
    		thisLayerNum --;
    		if(cur.left != null){
    			nextLayerNum ++;
    			q.offer(cur.left);
    		}
    		if(cur.right != null){
    			nextLayerNum ++;
    			q.offer(cur.right);
    		}
    		if(thisLayerNum == 0){
    			thisLayerNum = nextLayerNum;
    			nextLayerNum = 0;
    			result.add(thisLayerInt);
    			thisLayerInt = new ArrayList<Integer>();
    		}
    	}
    	return result;
    }
    
    //binary tree maximum path sum
    public int maxPathSum(TreeNode root) {
        // Start typing your Java solution below
        // DO NOT write main() function
    	int result = 0;
    	if(root == null)	return result;
    	combin1 res = findMaxSumPath(root);
    	result = res.maxSum;
    	return result;
    }
    class combin1{
    	int maxSum = 0;
    	int maxSumWithoutTurn = 0;
    }
    combin1 findMaxSumPath(TreeNode root){
    	combin1 result = new combin1();
    	if(root == null)	return result;
    	if(root.left == null && root.right == null){
    		result.maxSum = result.maxSumWithoutTurn = root.val;
    	}else if(root.right == null){
    		combin1 left = findMaxSumPath(root.left);
    		result.maxSumWithoutTurn = left.maxSumWithoutTurn > 0 ? left.maxSumWithoutTurn + root.val : root.val;
    		result.maxSum = result.maxSumWithoutTurn > left.maxSum ? result.maxSumWithoutTurn : left.maxSum;
    	}else if(root.left == null){
    		combin1 right = findMaxSumPath(root.right);
    		result.maxSumWithoutTurn = right.maxSumWithoutTurn > 0 ? right.maxSumWithoutTurn + root.val : root.val;
    		result.maxSum = result.maxSumWithoutTurn > right.maxSum ? result.maxSumWithoutTurn : right.maxSum;
    	}
    	else{
        	combin1 left = findMaxSumPath(root.left);
        	combin1 right = findMaxSumPath(root.right);
        	int tmp_max = left.maxSumWithoutTurn > right.maxSumWithoutTurn ? left.maxSumWithoutTurn + root.val : right.maxSumWithoutTurn + root.val;
        	result.maxSumWithoutTurn = tmp_max > root.val ? tmp_max : root.val;
        	int max_turn = left.maxSumWithoutTurn + right.maxSumWithoutTurn + root.val;
        	int local_max = result.maxSumWithoutTurn > max_turn ? result.maxSumWithoutTurn : max_turn;
        	int left_max = left.maxSum;
        	int right_max = right.maxSum;
        	if(left_max > right_max)	result.maxSum = local_max > left_max ? local_max : left_max;
        	else						result.maxSum = local_max > right_max ? local_max : right_max;
    	}
    	return result;
    }
    
    //binary tree zigzag level order traversal
    public ArrayList<ArrayList<Integer>> zigzagLevelOrder(TreeNode root) {
        // Start typing your Java solution below
        // DO NOT write main() function
    	ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
    	boolean left2right = true;
    	if(root == null)	return result;
    	Queue<TreeNode> q = new LinkedList<TreeNode>();
    	Stack<Integer> st = new Stack<Integer>();
    	ArrayList<Integer> al = new ArrayList<Integer>();
    	int numThisLayer = 1, numNextLayer = 0;
    	q.offer(root);
    	while(!q.isEmpty()){
    		TreeNode cur = q.poll();
    		numThisLayer --;
    		if(left2right)	al.add(cur.val);
    		else st.push(cur.val);
    		if(cur.left != null){
    			q.offer(cur.left);
    			numNextLayer ++;
    		}
    		if(cur.right != null){
    			q.offer(cur.right);
    			numNextLayer ++;
    		}
    		if(numThisLayer == 0){
    			numThisLayer = numNextLayer ;
    			numNextLayer = 0;
    			while(!left2right && !st.isEmpty())
    				al.add(st.pop());
    			result.add(al);
    			al = new ArrayList<Integer>();
    			left2right = left2right == true ? false : true;
    		}
    	}
    	return result;
    }

    //climing stairs 
    class matrix{
    	int a, b, c, d;
    	matrix(int aa, int bb, int cc, int dd) { a = aa; b = bb; c = cc; d = dd; }
    }
    matrix matrix_mul(matrix left, matrix right){
    	matrix res = new matrix(left.a * right.a + left.b * right.c, left.a * right.b + left.b * right.d, 
    							left.c * right.a + left.d * right.c, left.c * right.b + left.d * right.d);
    	return res;
    }
    matrix calFabonaci(int n, matrix[] record){
    	if(n == 1)	return record[1] = new matrix(1,1,1,0);
    	if(record[n] != null)	return record[n];
    	if(n%2 == 1){
    		matrix half = calFabonaci((n-1)/2,record);
    		return record[n] = matrix_mul(matrix_mul(half,half),new matrix(1,1,1,0));
    	}
    	else{
    		matrix half = calFabonaci(n/2,record);
    		return record[n] = matrix_mul(half,half);
    	}
    }
    public int climbStairs(int n) {
    	//logn Fabonaci
    	int result = 0;
    	if(n == 0 || n == 1 || n == 2) return n;
    	matrix[] record = new matrix[n+1];
    	matrix tmp_res = calFabonaci(n-2,record);
    	matrix tmp_res2 = matrix_mul(new matrix(2,1,1,0),tmp_res);
    	result = tmp_res2.a;
    	return result;
    }
    
    //Combination sum
    public void combinationSumRec(int[] candidates, int target, int k, ArrayList<Integer> cur, ArrayList<ArrayList<Integer>> result) {
    	if(target == 0){
    		ArrayList<Integer> tmp_res = new ArrayList<Integer>();
    		for(int i=0; i<cur.size(); i++)
    			tmp_res.add(cur.get(i));
    		result.add(tmp_res);
    		return;
    	}
    	if(target < 0 || k >= candidates.length)	return;
    	for(int i=0; i*candidates[k] <= target; i++){
    		if(i != 0)	cur.add(candidates[k]);
    		combinationSumRec(candidates, target-i*candidates[k], k+1, cur, result);
    	}
    	for(int i=1; i*candidates[k] <= target; i++)
    		cur.remove(cur.size()-1);
    }
    public ArrayList<ArrayList<Integer>> combinationSum(int[] candidates, int target) {
    	ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
    	ArrayList<Integer> cur = new ArrayList<Integer>();
    	Arrays.sort(candidates);
    	combinationSumRec(candidates, target, 0, cur, result);
    	return result;
    }

    //Combination sum2
    public void combinSum2Rec(int []num, int target, int k, ArrayList<Integer> cur, ArrayList<ArrayList<Integer>> result){
    	if(target == 0){
    		ArrayList<Integer> tmp_res = new ArrayList<Integer>();
    		for(int i=0;i<cur.size();i++)
    			tmp_res.add(cur.get(i));
    		result.add(tmp_res);
    		return;
    	}
    	if(target < 0 || k >= num.length)	return;
    	cur.add(num[k]);
    	combinSum2Rec(num,target-num[k],k+1,cur,result);
    	cur.remove(cur.size()-1);
    	int cur_val = num[k];
    	while(k<num.length && num[k] == cur_val)  k++;
    	combinSum2Rec(num,target,k,cur,result);
    }
    public ArrayList<ArrayList<Integer>> combinationSum2(int[] num, int target) {
    	ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
    	Arrays.sort(num);
    	ArrayList<Integer> cur = new ArrayList<Integer>();
    	combinSum2Rec(num,target,0,cur,result);
    	return result;
    }
    
    //combinations
    public void combinRec(int n, int k, int m, ArrayList<Integer> cur, ArrayList<ArrayList<Integer>> result){
    	if(k == 0){
    		ArrayList<Integer> tmp = new ArrayList<Integer>();
    		for(int i=0; i<cur.size(); i++)		tmp.add(cur.get(i));
    		result.add(tmp);
    		return;
    	}
    	for(int i=m;i<=n-k+1;i++){
    		cur.add(i);
    		combinRec(n,k-1,i+1,cur,result);
    		cur.remove(cur.size()-1);
    	}
    }
    public ArrayList<ArrayList<Integer>> combine(int n, int k) {
    	ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
    	ArrayList<Integer> cur = new ArrayList<Integer>();
    	combinRec(n,k,1,cur,result);
    	return result;
    }
    
    //construct binary tree from in order and post order traversal
    public TreeNode buildTreeRec(int []inorder, int []postorder, int s1, int e1, int s2, int e2){
    	TreeNode root = null;
    	if(s1>e1 || s2>e2)	return root;
    	root = new TreeNode(postorder[e2]);
    	int i = s1;
    	for(i=s1; i<=e1; i++){
    		if(inorder[i] == root.val)	break;
    	}
    	root.left = buildTreeRec(inorder, postorder, s1, i-1, s2, s2+(i-s1-1));
    	root.right = buildTreeRec(inorder, postorder, i+1, e1, s2+i-s1, e2-1);
    	return root;
    }
    public TreeNode buildTree(int[] inorder, int[] postorder) {
    	TreeNode root = null;
    	if(inorder.length == 0)	return root;
    	root = buildTreeRec(inorder, postorder, 0,inorder.length-1, 0,postorder.length-1);
    	return root;
    }

    //construct binary tree from pre order and in order traversal
    public TreeNode buildTree2Rec(int []preorder, int []inorder, int s1, int e1, int s2, int e2){
    	TreeNode root = null;
    	if(s1 > e1 || s2 > e2)	return root;
    	root = new TreeNode(preorder[s1]);
    	int i = s2;
    	for(i=s2; i<=e2; i++){
    		if(inorder[i] == root.val)	break;
    	}
    	root.left = buildTree2Rec(preorder, inorder, s1+1,s1+i-s2, s2,i-1);
    	root.right = buildTree2Rec(preorder, inorder, s1+i-s2+1,e1, i+1, e2);
    	return root;
    }
    public TreeNode buildTree2(int[] preorder, int[] inorder) {
        TreeNode root = null;
        if(preorder.length == 0)	return root;
        root = buildTree2Rec(preorder, inorder, 0,preorder.length-1, 0,inorder.length-1);
        return root;
    }
    
    //container with most water
    public int maxArea(int[] height) {
        int left = 0; 
        int right = height.length-1;
        int result = height[left] < height[right] ? height[left] * (height.length-1) : height[right] * (height.length-1);
        while(left < right){
        	int tmp_result = height[left] < height[right] ? height[left] * (right - left) : height[right] * (right - left) ;
        	result = tmp_result > result ? tmp_result : result;
        	if(height[left] < height[right]){
        		left ++;
        	}else{
        		right --;
        	}
        }	
        return result;
    }

    //convert sorted array to binary search tree
    public TreeNode array2bst(int []num, int s, int e){
    	TreeNode root = null;
    	if(s>e)	return root;
    	root = new TreeNode(num[(s+e)/2]);
    	root.left = array2bst(num,s,(s+e)/2-1);
    	root.right = array2bst(num,(s+e)/2+1,e);
    	return root;
    }
    public TreeNode sortedArrayToBST(int[] num) {
    	return array2bst(num,0,num.length-1);
    }
   
    //count and say
    public String countAndSay(int n) {
    	String result = "";
    	if(n==0)	return result;
    	result = "1";
    	String next;
    	for(int i=1;i<n;i++){
    		char []tmp = result.toCharArray();
    		next = "";
    		char prev = tmp[0];
    		int times = 1;
    		for(int j=1;j<tmp.length;j++){
    			if(tmp[j] == prev)	times ++;
    			else{
    				next += Integer.toString(times) + prev;
    				times = 1;
    				prev = tmp[j];
    			}
    		}
    		next += Integer.toString(times) + prev;
    		result = next;
    	}
    	return result;
    }

    //decode ways
    public int numDecodeRec(String s, int[]record, int k){
    	int result = 0;
    	if(k > s.length())	return result;
    	if(k == s.length())	return 1;
    	if(record[k] != -1)	return record[k];
    	if(k == s.length()-1){
    		return record[k] = (s.charAt(k) == '0' ? 0 : 1);
    	}
    	else{
    		char first = s.charAt(k);
    		char second = s.charAt(k+1);
    		if(first == '0')	return record[k] = 0;
    		if(first == '1')	return record[k] = (second == '0' ? numDecodeRec(s, record, k+2) : numDecodeRec(s,record,k+1) + numDecodeRec(s,record,k+2));
    		if(first == '2'){
    			if(second == '0')	return record[k] = numDecodeRec(s,record,k+2);
    			else if(second > '6')	return record[k] = numDecodeRec(s,record,k+1);
    			else	return record[k] = numDecodeRec(s,record,k+1) + numDecodeRec(s,record,k+2);
    		}
    		else return record[k] = (second == '0' ? 0 : numDecodeRec(s, record, k+1));
    	}
    }
    public int numDecodings(String s) {
        int result = 0;
        if(s == null || s.length() == 0)	return result;
        if(s.length() == 1)		return s.charAt(0) == '0' ? 0 : 1;
        int []record = new int[s.length()];
        for(int i=0; i<s.length(); i++)		record[i] = -1;
        result = numDecodeRec(s, record, 0);
        return result;
    }
    
    //distinct subsequences
    public int numDistinct(String S, String T) {
        int result = 0;
        if(S == null || T == null || S.length() == 0 || T.length() == 0)	return result;
        int [][]record = new int[S.length()+1][T.length()+1];
        for(int i=0;i<=T.length();i++)	record[S.length()][i] = 0;
        for(int i=0;i<=S.length();i++)	record[i][T.length()] = 1;
        for(int m=S.length()-1; m>=0;m--){
        	for(int n=T.length()-1; n>=0; n--){
        		if(S.charAt(m) == T.charAt(n))
        			record[m][n] = record[m+1][n+1] + record[m+1][n];
        		else
        			record[m][n] = record[m+1][n];
        	}
        }
        //result = numDistinctRec(S,T,record,0,0);
        result = record[0][0];
        return result;
    }
    
    //divide two integers
    public int divide(int dividend, int divisor) {
    	//search for a range divisor * 2^i < dividend < divisor * 2^(i+1) first
    	//then binary search (record 1,2,4...2^i)
    	//to avoid overflow, turn positive int to negative 
    	boolean isNeg = false;
    	if(dividend > 0){
    		if(divisor < 0)
    			isNeg = true;
    		else
    			divisor = -divisor;
    		dividend = -dividend;
    	}else if(divisor > 0){
    		isNeg = true;
    		divisor = -divisor;
    	}
    	int []record = new int[66];
    	int index = 0, sum = divisor;
    	while(sum > dividend){
    		record[index ++] = sum; 
    		sum += sum;
    		if(sum >= 0)	break;
    	}
    	if(sum == dividend)	return  isNeg ? -(1<<index) : (1<<index);
    	if(index == 0)	return 0;
    	int result = (1 << (index-1)), rest = dividend - record[index-1];
    	index --;
    	while(rest < 0 && index > 0){
    		if(rest <= record[index-1]){
    			rest -= record[index-1];
    			result += (1 << (index - 1));
    		}
    		index --;
    	}
    	if(isNeg)
    		return -result;
    	else
    		return result;
    }

    //edit distance
    public int minDistance(String word1, String word2) {
    	int result = 0;
    	int [][]record = new int[word1.length()+1][word2.length()+1];
    	for(int i=1; i<=word1.length(); i++)	record[i][0] = i;
    	for(int i=1; i<=word2.length(); i++)	record[0][i] = i;
    	record[0][0] = 0;
    	for(int i=1; i<=word1.length(); i++){
    		for(int j = 1; j<=word2.length(); j++){
    			int min1 = record[i-1][j] < record[i][j-1] ? record[i-1][j]+1 : record[i][j-1]+1;
    			if(word1.charAt(i-1) == word2.charAt(j-1)){
    				int min2 = min1 < record[i-1][j-1] ? min1 : record[i-1][j-1];
    				record[i][j] = min2;
    			}else{
    				int min2 = min1 < record[i-1][j-1] + 1 ?  min1 : record[i-1][j-1]+1;
    				record[i][j] = min2;
    			}
    		}
    	}
    	result = record[word1.length()][word2.length()];
    	return result;
    }

    //first missing positive
    public int firstMissingPositive(int[] A) {
        int result = 0;
        for(int i=0; i<A.length; i++){
        	if(A[i] <= 0 || A[i] > A.length)	A[i] = 0;
        	else{
        		int prev = A[i];
        		A[i] = 0;
        		while(true){
        			if(A[prev-1] <= 0 || A[prev-1] > A.length){
        				A[prev-1] = prev;
        				break;
        			}else if(A[prev-1] != prev){
        				int tmp = A[prev-1];
        				A[prev-1] = prev;
        				prev = tmp;
        			}else break;
        		}
        	}
        }
        int i = 0;
        for(i=0; i<A.length; i++){
        	if(A[i] == 0)	break;
        }
        result = i+1;
        return result;
    }

    //flatten binary tree to linked list
    public void flatten(TreeNode root) {
    	if(root == null)	return;
    	Stack<TreeNode> st = new Stack<TreeNode>();
    	TreeNode left = root.left, right = root.right;
    	while(left != null || right != null || !st.isEmpty()){
    		if(left == null && right == null){
    			TreeNode next = st.pop();
    			root.right = next;
    			root = next;
    		}else{
    			if(right == null){
    				root.right = left;
    				root.left = null;
    				root = left;
    			}else{
    				if(left == null)
    					root = right;
    				else{
    					st.push(right);
    					root.right = left;
    					root.left = null;
    					root = left;
    				}
    			}
    		}
    		left = root.left;
    		right = root.right;
    	}
    }

    //generate parenthesis
    public void generatePathRec(int left,int right, char [] cur, int n, ArrayList<String> result){
    	if(right == n){
    		result.add(new String(cur));
    		return;
    	}
    	if(left == right){
    		cur[left+right] = '(';
    		generatePathRec(left+1,right,cur,n,result);
    	}
    	else{
    		if(left < n){
    			cur[left+right] = '(';
        		generatePathRec(left+1,right,cur,n,result);
    		}
			cur[left+right] = ')';
			generatePathRec(left,right+1,cur,n,result);
    	}
    }
    public ArrayList<String> generateParenthesis(int n) {
    	ArrayList<String> result = new ArrayList<String>();
    	char []cur = new char[n*2];
    	generatePathRec(0,0,cur,n,result);
    	return result;
    }
    
    //gray code
    public ArrayList<Integer> grayCode(int n) {
        ArrayList<Integer> result = new ArrayList<Integer>();
        int []record = new int[1<<n];
        record[0] = 0; 
        for(int i=1;i<=n;i++){
        	int cnt = 1<<(i-1);
        	for(int j=(1<<(i-1))-1; j>=0 ;j--)
        		record[cnt++] = (1<<(i-1)) | record[j];
        }
        for(int i=0;i<(1<<n);i++)
        	result.add(record[i]);
        return result;
    }

    //implement strstr()
    //need to implement KMP later
    public String strStr(String haystack, String needle) {
        String result = null;
        for(int i=0; i<=haystack.length()-needle.length(); i++){
        	int j = 0;
        	for(j=0; j<needle.length(); j++){
        		if(haystack.charAt(i+j) != needle.charAt(j))
        			break;
        	}
        	if(j==needle.length())
        		return haystack.substring(i);
        }
        return result;
    }
    
    //insert interval
    public ArrayList<Interval> insert(ArrayList<Interval> intervals, Interval newInterval) {
        ArrayList<Interval> result = new ArrayList<Interval>();
        if(intervals.size() == 0){
            result.add(newInterval);
        	return result;
        }
        boolean isNewintervalInserted = false;
        Interval replaceInterval = new Interval(newInterval.start, newInterval.end);
        for(int i=0; i<intervals.size(); i++){
        	Interval cur = intervals.get(i);
        	if(replaceInterval.start >= cur.start && replaceInterval.start <= cur.end){
        		if(replaceInterval.end <= cur.end){
        			isNewintervalInserted = true;
        			result.add(cur);
        		}else{
        			replaceInterval.start = cur.start;
        		}
        	}else if(replaceInterval.end >= cur.start && replaceInterval.end <= cur.end){
        		if(replaceInterval.start >= cur.start){
        			isNewintervalInserted = true;
        			result.add(cur);
        		}else{
        			replaceInterval.end = cur.end;
        			isNewintervalInserted = true;
        			result.add(replaceInterval);
        		}
        	}else if(replaceInterval.end < cur.start || replaceInterval.start > cur.end){
        		if(!isNewintervalInserted && cur.start > replaceInterval.end){
        			isNewintervalInserted = true;
        			result.add(replaceInterval);
        		}
        		result.add(cur);
        	}
        }
        if(!isNewintervalInserted)	result.add(replaceInterval);
        return result;
    }

    //interleaving string
    public boolean isInterleave(String s1, String s2, String s3) {
        // Start typing your Java solution below
        // DO NOT write main() function
        boolean result = false;
        if(s1 == null || s2 == null || s3 == null)	return result;
        if(s3.length() != s1.length() + s2.length())	return false;
        boolean [][]record = new boolean[s1.length()+1][s2.length()+1];
        record[0][0] = true;
        for(int i=1; i<=s1.length(); i++)	record[i][0] = record[i-1][0] && (s1.charAt(i-1) == s3.charAt(i-1));
        for(int i=1; i<=s2.length(); i++)	record[0][i] = record[0][i-1] && (s2.charAt(i-1) == s3.charAt(i-1));
        for(int i=1; i<=s1.length(); i++){
        	for(int j=1; j<=s2.length(); j++){
        		record[i][j] = (record[i-1][j] && (s1.charAt(i-1) == s3.charAt(i+j-1)) || (record[i][j-1] &&(s2.charAt(j-1) == s3.charAt(i+j-1))));
        	}
        }
        result = record[s1.length()][s2.length()];
        return result;
    }

    //jump game
    public boolean canJump(int[] A) {
        if(A.length == 0)	return true;
        int cur_max = A[0], i=1;
        for(i=1; i<A.length && i <= cur_max; i++){
        	if(i+A[i] > cur_max)
        		cur_max = i + A[i];
        }
        if(i == A.length)	return true;
        return false;
    }

    //jump game 2
    public int jump(int[] A) {
        int result = 0;
        if(A.length == 0)	return 0;
        int []record = new int[A.length];
        record[0] = 0;
        int cur_max = 0;
        for(int i=0;i<A.length; i++){
        	if(i+A[i] > cur_max){
        		for(int j=cur_max+1;j<=i+A[i] && j < record.length ;j++)
        			record[j] = record[i] + 1;
        		cur_max = i+A[i];
        	}
        }
        result = record[A.length-1];
        return result;
    }

    //largest rectangle in histogram
    public int largestRectangleArea(int[] height) {
        int result = 0;
        if(height.length == 0)	return 0;
        int []leftX = new int[height.length], leftHeight = new int[height.length];
        int []rightX = new int[height.length], rightHeight = new int[height.length];
        leftX[0] = -1; leftHeight[0] = height[0]; rightX[height.length-1] = height.length; rightHeight[height.length-1] = height[height.length-1];
        for(int i=1; i<height.length; i++){
        	if(height[i-1] < height[i]){
        		leftX[i] = i-1;
        		leftHeight[i] = height[i-1];
        	}else{
        		int tmpHeight = leftHeight[i-1];
        		int tmpX = leftX[i-1];
        		while(tmpHeight >= height[i]){
        			if(tmpX == -1)	{
        				tmpHeight = height[i];
        				break;
        			}
        			tmpHeight = leftHeight[tmpX];
        			tmpX = leftX[tmpX];
        		}
        		leftX[i] = tmpX;
        		leftHeight[i] = tmpHeight;
        	}
        }
        for(int i=height.length-2; i>=0 ;i--){
        	if(height[i+1] < height[i]){
        		rightX[i] = i+1;
        		rightHeight[i] = height[i+1];
        	}else{
        		int tmpHeight = rightHeight[i+1];
        		int tmpX = rightX[i+1];
        		while(tmpHeight >= height[i]){
        			if(tmpX == height.length){
        				tmpHeight= height[i];
        				break;
        			}
        			tmpHeight = rightHeight[tmpX];
        			tmpX = rightX[tmpX];
        		}
        		rightX[i] = tmpX;
        		rightHeight[i] = tmpHeight;
        	}
        }
        for(int i=0; i<height.length; i++){
        	if(height[i] * (rightX[i] - leftX[i] - 1) > result)
        		result = height[i] * (rightX[i] - leftX[i] - 1);
        }
        return result;
    }
    
    //length of last word
    public int lengthOfLastWord(String s) {
    	// "  a b   "
        int result = 0;
        if(s == null || s.length() == 0)	return 0;
        boolean startCnt = false;
        int cnt = 0 ;
        for(int i=s.length()-1;i>=0;i--){
        	if(s.charAt(i) == ' '){
        		if(!startCnt && cnt == 0)	startCnt = true;
        		else if(cnt != 0)	break;
        	}
        	else{
        		cnt ++;
        	}
        }
        result = cnt;
        return result;
    }

    //letter combinations of phone number
    public void letterRec(String digits, int cur,char []tmp, ArrayList<String> result){
    	if(cur == digits.length()){
    		result.add(new String(tmp));
    		return;
    	}
    	if(digits.charAt(cur) >= '2' && digits.charAt(cur) <= '6'){
    		int d = digits.charAt(cur) - '2';
    		char c1 = (char)('a' + d*3);
    		char c2 = (char)('a' + d*3+1);
    		char c3 = (char)('a' + d*3+2);
    		tmp[cur] = c1;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c2;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c3;
    		letterRec(digits,cur+1,tmp,result);
    	}else if(digits.charAt(cur) == '7'){
    		char c1 = 'p', c2 = 'q', c3 = 'r', c4 = 's';
    		tmp[cur] = c1;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c2;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c3;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c4;
    		letterRec(digits,cur+1,tmp,result);
    	}else if(digits.charAt(cur) == '8'){
    		char c1 = 't', c2 = 'u', c3 = 'v';
    		tmp[cur] = c1;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c2;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c3;
    		letterRec(digits,cur+1,tmp,result);
    	}else if(digits.charAt(cur) == '9'){
    		char c1 = 'w', c2 = 'x', c3 = 'y', c4 = 'z';
    		tmp[cur] = c1;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c2;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c3;
    		letterRec(digits,cur+1,tmp,result);
    		tmp[cur] = c4;
    		letterRec(digits,cur+1,tmp,result);
    	}else if(digits.charAt(cur) == '0'){
    		tmp[cur] = ' ';
    		letterRec(digits,cur+1,tmp,result);
    	}
    }
    public ArrayList<String> letterCombinations(String digits) {
        ArrayList<String> result = new ArrayList<String>();
        if(digits == null || digits.length() == 0)	{
        	result.add("");
        	return result;
        }
        char []tmp = new char[digits.length()];
        letterRec(digits, 0, tmp, result);
        return result;
    }
    
    //longest common prefix
    public String longestCommonPrefix(String[] strs) {
        String result = "";
        if(strs.length == 0 || strs[0].length() == 0)	return result;
        for(int index = 0; index < strs[0].length(); index ++){
        	char c = strs[0].charAt(index);
        	int i = 0;
        	for(i=0;i<strs.length; i++){
        		if(index >= strs[i].length() || strs[i].charAt(index) != c)
        			break;
        	}
        	if(i == strs.length){
        		result = result + c;
        	}else	break;
        }
        return result;
    }

    //longest consecutive sequence
    public int longestConsecutive(int[] num2) {
        int result = 1;
        if(num2 == null || num2.length == 0)	return result;
        int []num = new int[num2.length];
        HashSet<Integer> hs = new HashSet<Integer>();
        int cnt = 0;
        for(int i=0; i<num2.length; i++){
        	if(!hs.contains(num2[i])){
        		hs.add(num2[i]);
        		num[cnt++] = num2[i];
        	}
        }
        HashMap<Integer, Integer> upperBound = new HashMap<Integer,Integer>();
        HashMap<Integer, Integer> lowerBound = new HashMap<Integer,Integer>();
        for(int i=0; i<cnt; i++){
        	if(!lowerBound.containsKey(num[i]+1) && !upperBound.containsKey(num[i]-1)){
        		upperBound.put(num[i], 1);
        		lowerBound.put(num[i], 1);
        	}else{
        		if(lowerBound.containsKey(num[i]+1) && !upperBound.containsKey(num[i]-1)){
        			int oriLength = lowerBound.get(num[i]+1);
        			result = oriLength + 1 > result ? oriLength + 1 : result;
        			lowerBound.put(num[i], oriLength + 1);
        			upperBound.put(num[i] + oriLength, oriLength + 1);
        			lowerBound.remove(num[i]+1);
        		}
        		else if(!lowerBound.containsKey(num[i]+1) && upperBound.containsKey(num[i]-1)){
        			int oriLength = upperBound.get(num[i]-1); 
        			result = oriLength + 1 > result ? oriLength + 1 : result;
        			upperBound.put(num[i], oriLength + 1);
        			lowerBound.put(num[i]-oriLength, oriLength+1);
        			upperBound.remove(num[i]-1);
        		}
        		else{
        			int oriLenLow = lowerBound.get(num[i]+1);
        			int oriLenHigh = upperBound.get(num[i]-1);
        			result = oriLenLow + oriLenHigh + 1 > result ? oriLenLow + oriLenHigh + 1 : result;
        			lowerBound.put(num[i] - oriLenHigh, oriLenLow+oriLenHigh+1);
        			upperBound.put(num[i] + oriLenLow, oriLenLow+oriLenHigh+1);
        			lowerBound.remove(num[i]+1);
        			upperBound.remove(num[i]-1);
        		}
        	}
        }
        return result;
    }
        
    //longest palindrome substring
    //need to know O(N) solution later
    public String longestPalindrome(String s) {
    	String result = "";
    	if(s== null || s.length() == 0)	return result;
    	int max_len = 0;
    	for(int i=0; i<s.length(); i++){
    		char center = s.charAt(i);
    		int tmp_len = 1;
    		int j = 1;
    		for(j=1; i-j>=0 && i+j<s.length();j++){
    			if(s.charAt(i-j) == s.charAt(i+j)){
    				tmp_len += 2; 
    			}else{
    				break;
    			}
    		}
    		if(tmp_len > max_len){
    			max_len = tmp_len;
    			result = s.substring(i-j+1,i+j);
    		}
    		if(i > 0 && s.charAt(i-1) == center){
    			j = 1; tmp_len = 2;
    			for(j=1; i-1-j >= 0 && i+j < s.length(); j++){
    				if(s.charAt(i-1-j) == s.charAt(i+j))
    					tmp_len += 2;
    				else
    					break;
    			}
    			if(tmp_len > max_len){
    				max_len = tmp_len;
    				result = s.substring(i-j,i+j);
    			}
    		}
    	}
    	return result;
    }
    
    //longest substring without repeating characters
    public int lengthOfLongestSubstring(String s) {
        int result = 0;
        if(s == null || s.length() == 0)	return 0;
        HashSet<Character> hs = new HashSet<Character>();
        int start = 0, end = 0;
        for(int i = 0; i<s.length(); i++){
        	if(! hs.contains(s.charAt(i))){
        		hs.add(s.charAt(i));
        		end = i;
        		result = end-start+1 > result ? end-start+1 : result;
        	}else{
        		int j=start;
        		while(s.charAt(j) != s.charAt(i))
        			hs.remove(s.charAt(j++));
        		j ++;
        		start = j;
        		end = i;
        	}
        }
        return result;
    }

    //longest valid parentheses
    //need a simple solution
    public int longestValidParentheses(String s) {
    	//need to look back for 2 steps 
        int result = 0;
        if(s == null || s.length() == 0 || s.length() == 1)	return 0;
        class tmp_class{
        	char c;
        	int x;
        	tmp_class(char cc, int xx) { c = cc; x = xx; }
        }
        Stack<tmp_class> st = new Stack<tmp_class>();
        for(int i=0; i<s.length(); i++){
        	char c = s.charAt(i);
        	if(c == '('){
        		st.push(new tmp_class('(',i));
        	}else{
        		if(!st.isEmpty()){
        			if(st.peek().c == '('){
        				st.pop();
        				if(!st.isEmpty() && st.peek().c == 'a'){
        					result = i-st.peek().x + 1 > result ? i-st.peek().x + 1 : result;
        				}else{
        					st.push(new tmp_class('a',i-1));
        					result = result > 2 ? result : 2;
        				}
        			}else{
        				st.pop();
        				if(!st.isEmpty() && st.peek().c == '('){
        					int prev_i = st.pop().x;
        					if(st.isEmpty() || st.peek().c == '('){
        						st.push(new tmp_class('a',prev_i));
        						result = i-prev_i+1 > result ? i-prev_i+1:result;
        					}else{
        						int prev_i2 = st.pop().x;
        						st.push(new tmp_class('a',prev_i2));
        						result = i-prev_i2+1 > result ? i-prev_i2+1:result;
        					}
        				}
        				else{
        					st.clear();
        				}
        			}
        		}
        	}
        }
        return result;
    }

    //maximal rectangle
    //rewrite
    public int maximalRectangle(char[][] matrix) {
        int result = 0;
        if(matrix.length == 0 || matrix[0].length == 0)
        	return 0;
        int []prev = new int[matrix[0].length];
        int []height = new int[matrix[0].length];
        
        for(int i=0; i<matrix.length; i++){
        	for(int j=0; j<matrix[0].length; j++){
        		if(matrix[i][j] == '1')
        			height[j] = prev[j] + 1;
        		else
        			height[j] = 0;
        	}
        	
        	result = largestRectangleArea(height) > result ? largestRectangleArea(height) : result; 
        	
        	for(int j=0; j<matrix[0].length; j++)
        		prev[j] = height[j];
        }
        return result; 
    }
    
    //maximum depth of binary tree
    public int maxDepth(TreeNode root) {
        if(root == null)	return 0;
        int left = maxDepth(root.left);
        int right = maxDepth(root.right);
        return left > right ? left + 1 : right + 1;
    }
    
    //maximum subarray
    public int maxSubArray(int[] A) {
        int result = Integer.MIN_VALUE;
        if(A.length == 0)	return 0;
        int cur_sum = 0;
        for(int i=0; i<A.length; i++){
        	if(cur_sum < 0){
        		cur_sum = A[i];
        		result = A[i] > result ? A[i] : result;
        	}else{
        		cur_sum += A[i];
        		result = cur_sum > result ? cur_sum : result;
        	}
        }
        return result;
    }

    //merge intervals
    public ArrayList<Interval> merge(ArrayList<Interval> intervals) {
        ArrayList<Interval> result = new ArrayList<Interval>();
        if(intervals.size() == 0)	return result;
        if(intervals.size() == 1)	{ result.add(intervals.get(0)); return result; }
        Collections.sort(intervals, new Comparator<Interval>(){
			@Override
			public int compare(Interval o1, Interval o2) {
				return o1.start - o2.start;
			}
        });
        Interval prev = intervals.get(0); 
        for(int i=1; i<intervals.size(); i++){
        	Interval cur = intervals.get(i);
        	if(cur.start > prev.end){
        		result.add(prev);
        		prev = cur;
        	}else{
        		if(cur.end > prev.end)
        			prev.end = cur.end;
        	}
        }
        result.add(prev);
        return result;
    }

    //merge k sorted lists
    public ListNode mergeKLists(ArrayList<ListNode> lists) {
        ListNode result = null, cur = null;
        if(lists.size() == 0)	return result;
        PriorityQueue<ListNode> min_heap = new PriorityQueue<ListNode>(lists.size(),new Comparator<ListNode>(){
			@Override
			public int compare(ListNode o1, ListNode o2) {
				return o1.val-o2.val;
			}
        });
        for(int i=0; i<lists.size(); i++){
        	if(lists.get(i) != null)
        		min_heap.offer(lists.get(i));
        }
        while(!min_heap.isEmpty()){
        	ListNode cur_min = min_heap.poll();
        	if(result == null){
        		result = cur = cur_min;
        	}else{
        		cur.next = cur_min;
        		cur = cur_min;
        	}
        	if(cur_min.next != null)
        		min_heap.offer(cur_min.next);
        }
        return result;
    }

    //merge sorted lists
    public void merge(int A[], int m, int B[], int n) {
        if(n == 0)	return;
        if(m == 0) {
        	for(int i=0;i<n;i++)
        		A[i] = B[i];
        	return;
        }
    	int a = m-1, b = n-1;
        while(a>=0 && b>=0){
        	if(A[a] >= B[b]){
        		A[a+b+1] = A[a];
        		a --;
        	}else{
        		A[a+b+1] = B[b];
        		b--;
        	}
        }
        for(;b>=0;b--)
        	A[b] = B[b];
    }

    //merge two sorted lists
    public ListNode mergeTwoLists(ListNode l1, ListNode l2) {
        ListNode res = null, cur = null;
        if(l1 == null && l2 == null)	return null;
        if(l1 == null)	return l2;
        if(l2 == null)	return l1;
        while(l1 != null && l2 != null){
        	ListNode tmp = null;
        	if(l1.val < l2.val){
        		tmp = l1;
        		l1 = l1.next;
        	}else{
        		tmp = l2;
        		l2 = l2.next;
        	}
        	if(res == null)
        		res = cur = tmp;
        	else{
        		cur.next = tmp;
        		cur = tmp;
        	}
        }
        while(l1 != null){
        	cur.next = l1;
        	cur = l1;
        	l1 = l1.next;
        }
        while(l2 != null){
        	cur.next = l2;
        	cur = l2;
        	l2 = l2.next;
        }
        return res;
    }

    //minimum depth of binary tree
    public int minDepth(TreeNode root) {
        if(root == null)	return 0;
        if(root.left == null && root.right == null)	return 1;
        if(root.left == null)	return minDepth(root.right) + 1;
        if(root.right == null)	return minDepth(root.left) + 1;
        int tmpl, tmpr;
        return (tmpl = minDepth(root.left)) < (tmpr = minDepth(root.right)) ? tmpl + 1 : tmpr + 1;   
    }

    //minimum path sum
    public int minPathSum(int[][] grid) {
    	if(grid.length == 0 || grid[0].length == 0)	return 0;
        int []prev = new int[grid[0].length], cur;
        prev[0] = grid[0][0];
        for(int i=1; i<grid[0].length; i++)
        	prev[i] = grid[0][i] + prev[i-1];
        cur = new int[grid[0].length];
        for(int i=1; i<grid.length; i++){
        	cur[0] = prev[0] + grid[i][0];
        	for(int j=1; j<grid[0].length; j++)
        		cur[j] = prev[j] < cur[j-1] ? prev[j] + grid[i][j] : cur[j-1] + grid[i][j];
        	prev = cur;
        	cur = new int[grid[0].length];
        }
        return prev[grid[0].length-1];
    }

    //minimum window substring
    public String minWindow(String S, String T) {
    	//shirnk the left side every time when find a new matchs
        if(S.length() == 0 || T.length() == 0 || S.length() < T.length())	return "";
        String result = "";
        HashMap<Character, Integer> hm = new HashMap<Character, Integer>();
        HashMap<Character, Integer> dict = new HashMap<Character, Integer>();
        Queue<Integer> q = new LinkedList<Integer>();
        for(int i=0; i<T.length(); i++){
        	if(!dict.containsKey(T.charAt(i)))
        		dict.put(T.charAt(i), 1);
        	else
        		dict.put(T.charAt(i), dict.get(T.charAt(i))+1);
        }
        int total_num = 0;
        int min_len = S.length();
        for(int i=0; i<S.length(); i++){
        	char c = S.charAt(i);
        	if(dict.containsKey(c)){
        		q.offer(i);
        		if(!hm.containsKey(c)){
        			hm.put(c, 1);
        			total_num ++;
        		}
        		else{
        			hm.put(c, hm.get(c)+1);
        			if(hm.get(c) <= dict.get(c)){
        				total_num ++;
        			}
        		}
        		if(total_num == T.length()){
        			//shrink
        			while(!q.isEmpty()){
        				int tmp = q.peek();
        				char tmpc = S.charAt(tmp);
        				if(hm.get(tmpc) > dict.get(tmpc)){
        					q.poll();
        					hm.put(tmpc, hm.get(tmpc)-1);
        				}else{
        					if(hm.get(tmpc) == 1)
        						hm.remove(tmpc);
        					else
        						hm.put(tmpc, hm.get(tmpc)-1);
        					total_num --;
        					break;
        				}
        			}
        			if(i-q.peek()+1 <= min_len){
        				min_len = i-q.peek()+1;
        				result = S.substring(q.peek(),i+1);
        			}
        			q.poll();
        		}
        	}
        }
        return result;
    }
    
    //N-queens
    public void solveNQRec(int n, int cur, boolean []vert, boolean []left, boolean []right, ArrayList<String[]> result, String []res){
    	if(cur == n){
    		String []tmp = new String[n];
    		for(int i=0; i<res.length; i++)
    			tmp[i] = res[i];
    		result.add(tmp);
    		return;
    	}
    	for(int i=0;i<n;i++){
    		if(vert[i] == false && left[cur-i+n-1] == false && right[cur+i] == false){
    			vert[i] = left[cur-i+n-1] = right[cur+i] = true;
    			String tmpS = "";
    			for(int j=0; j<i; j++)	tmpS += '.';
    			tmpS += 'Q';
    			for(int j=i+1;j<n;j++)	tmpS += '.';
    			res[cur] = tmpS;
    			solveNQRec(n,cur+1,vert,left,right,result,res);
    			vert[i] = left[cur-i+n-1] = right[cur+i] = false;
    		}
    	}
    }
    public ArrayList<String[]> solveNQueens(int n) {
    	ArrayList<String []> result = new ArrayList<String[]>();
    	if(n == 0) return result;
        boolean []vert = new boolean[n];
        boolean []left = new boolean[2*n-1], right = new boolean[2*n-1];
        String []res = new String[n];
        solveNQRec(n,0,vert,left,right,result,res);
        return result;
    }

    //N-queens 2
    public int solveNQNumRec(int n, int cur, boolean []vert, boolean []left, boolean []right){
    	if(cur == n)	return 1;
    	int num = 0;
    	for(int i=0; i<n; i++){
    		if(vert[i] == false && left[cur-i+n-1] == false && right[cur+i] == false){
    			vert[i] = left[cur-i+n-1] = right[cur+i] = true;
    			num += solveNQNumRec(n,cur+1,vert,left,right);
    			vert[i] = left[cur-i+n-1] = right[cur+i] = false;
    		}
    	}
    	return num;
    }
    public int totalNQueens(int n) {
        if(n==0)	return 0;
        boolean []vert = new boolean[n];
        boolean []left = new boolean[2*n-1], right = new boolean[2*n-1];
        return solveNQNumRec(n,0,vert,left,right);
    }
    
    //next permutation
    public void nextPermutation(int[] num) {
        if(num.length == 1 || num.length == 0)	return;
        int i = num.length-1, start = num.length-1, end = num.length-1;
        for(; i>0;i--){
        	if(num[i] > num[i-1]){
        		break;
        	}
        }
        start = i;
        if(start != 0){
        	for(i=start;i<num.length;i++){
        		if(num[i] <= num[start-1])
        			break;
        	}
        	int tmp = num[i-1];
        	num[i-1] = num[start-1];
        	num[start-1] = tmp;
        }
        while(start<end){
        	int tmp = num[end];
        	num[end] = num[start];
        	num[start] = tmp;
        	start ++;
        	end --;
        }
    }

    //palindrome number
    public boolean isPalindrome(int x) {
        if(x<0)	return false;
        if(x==0)	return true;
        int bits = 0;
        int tmp = x;
        while(tmp != 0){
        	tmp /= 10;
        	bits ++;
        }
        int Highbase = 1;
        for(int i=1;i<bits;i++)
        	Highbase *= 10;
        while(x != 0){
        	if(x % 10 != x / Highbase)
        		return false;
        	x %= Highbase;
        	x /= 10;
        	Highbase /= 100;
        }
        return true;
    }
    
    //palindrome partition
    public void searchPartRec(String s, int n, boolean [][]record, ArrayList<String> tmp, ArrayList<ArrayList<String>> result){
    	if(n==s.length()){
    		ArrayList<String> tmp_res = new ArrayList<String>();
    		for(int i=0;i<tmp.size();i++)
    			tmp_res.add(tmp.get(i));
    		result.add(tmp_res);
    		return;
    	}
    	for(int i=n; i<s.length(); i++){
    		if(record[n][i]){
    			tmp.add(s.substring(n,i+1));
    			searchPartRec(s,i+1,record,tmp,result);
    			tmp.remove(tmp.size()-1);
    		}
    	}
    }
    public ArrayList<ArrayList<String>> partition(String s) {
        ArrayList<ArrayList<String>> result = new ArrayList<ArrayList<String>>();
        if(s == null || s.length() == 0)	return result;
        boolean [][]record = new boolean [s.length()][s.length()];
        for(int i=0;i<s.length();i++)
        	record[i][i] = true;
        for(int i=0; i<s.length()-1; i++)
        	record[i][i+1] = s.charAt(i) == s.charAt(i+1) ? true : false;
        for(int d=2;d<s.length();d++){
        	for(int i=0;i+d<s.length();i++){
        		record[i][i+d] = record[i+1][i+d-1] && (s.charAt(i) == s.charAt(i+d)) ? true : false;
        	}
        }
        ArrayList<String> tmp = new ArrayList<String>();
        searchPartRec(s,0,record,tmp,result);
        return result;
    }
    
    //palindrome partition 2
    //BFS of graph
    public int minCut(String s) {
        if(s == null || s.length() == 0)	return 0;
        boolean [][]record = new boolean [s.length()][s.length()];
        for(int i=0; i<s.length(); i++)
        	record[i][i] = true;
        for(int i=0; i<s.length()-1; i++)
        	record[i][i+1] = s.charAt(i) == s.charAt(i+1) ? true : false;
        for(int d=2; d<s.length(); d++){
        	for(int i=0; i+d<s.length(); i++){
        		record[i][i+d] = record[i+1][i+d-1] && (s.charAt(i) == s.charAt(i+d)) ? true : false;
        	}
        }
        //BFS of graph
        int step = 0, cur_num = 0, prev_num = 0;
        boolean []flag = new boolean[s.length()];
        Queue<Integer> q = new LinkedList<Integer>();
        q.offer(-1);
        prev_num = 1;
        cur_num = 0;
        while(!q.isEmpty()){
        	int cur_node = q.poll();
        	prev_num --;
        	if(cur_node == s.length()-1)	return step-1;
        	for(int i=cur_node+1;i<s.length();i++){
        		if(record[cur_node+1][i] == true && flag[i] == false){
        			q.offer(i);
        			flag[i] = true;
        			cur_num ++;
        		}
        	}
        	if(prev_num == 0){
        		step ++;
        		prev_num = cur_num; 
        		cur_num = 0;
        	}
        }
        return step;
    }
    
    //partition list
    public ListNode partition(ListNode head, int x) {
    	//mention the last pointer
        ListNode result = null;
        if(head == null)	return null;
        ListNode leftHead = null, left = null,  rightHead = null, right = null;
        while(head != null){
        	if(head.val < x){
        		if(leftHead == null)
        			leftHead = left = head;
        		else{
        			left.next = head;
        			left = head;
        		}
        	}else{
        		if(rightHead == null)
        			rightHead = right = head;
        		else{
        			right.next = head;
        			right = head;
        		}
        	}
        	head = head.next;
        }
        if(right != null)	right.next = null;
        if(left != null){
        	left.next = rightHead;
        	result = leftHead;
        }
        else
        	result = rightHead;
        return result;
    }

    //pascal's triangle
    public ArrayList<ArrayList<Integer>> generate(int numRows) {
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        if(numRows <= 0)	return result;
        ArrayList<Integer> row = new ArrayList<Integer>();
        row.add(1);
        result.add(row);
        for(int i=2;i<=numRows;i++){
        	ArrayList<Integer> prev = result.get(i-2);
        	ArrayList<Integer> cur = new ArrayList<Integer>();
        	for(int j=0;j<i;j++){
        		if(j==0 || j==i-1)
        			cur.add(1);
        		else
        			cur.add(prev.get(j-1)+prev.get(j));
        	}
        	result.add(cur);
        }
        return result;
    }
    
    //pascal's triangle 2
    public ArrayList<Integer> getRow(int rowIndex) {
        ArrayList<Integer> result = new ArrayList<Integer>();
        if(rowIndex < 0)	return result;
        int []prev = new int[rowIndex+1];
        int []cur = new int[rowIndex+1];
        prev[0] = 1;
        for(int i=2;i<=rowIndex+1;i++){
        	for(int j=0;j<i;j++){
        		if(j==0 || j==i-1)
        			cur[j] = 1;
        		else
        			cur[j] = prev[j-1]+prev[j];
        	}
        	for(int j=0;j<i;j++)
        		prev[j] = cur[j];
        }
        for(int i=0;i<rowIndex+1;i++)
        	result.add(prev[i]);
        return result;
    }

    //path sum
    public boolean hasPathSum(TreeNode root, int sum) {
    	if(root == null)	return false;
    	if(root.left == null && root.right == null)
    		return root.val == sum ? true : false;
    	if(root.left == null)
    		return hasPathSum(root.right,sum-root.val);
    	if(root.right == null)
    		return hasPathSum(root.left, sum-root.val);
    	if(hasPathSum(root.left,sum-root.val))
    		return true;
    	return hasPathSum(root.right,sum-root.val);
    }

    //path sum 2
    public void pathSumRec(TreeNode root, int sum, ArrayList<Integer> cur, ArrayList<ArrayList<Integer>> result){
    	if(root.left == null && root.right == null && sum == root.val){
    		cur.add(root.val);
    		ArrayList<Integer> tmp_res = new ArrayList<Integer>();
    		for(int i=0;i<cur.size();i++)
    			tmp_res.add(cur.get(i));
    		result.add(tmp_res);
    		cur.remove(cur.size()-1);
    		return;
    	}
    	if(root.left != null){
    		cur.add(root.val);
    		pathSumRec(root.left,sum-root.val,cur,result);
    		cur.remove(cur.size()-1);
    	}
    	if(root.right != null){
    		cur.add(root.val);
    		pathSumRec(root.right,sum-root.val,cur,result);
    		cur.remove(cur.size()-1);
    	}
    }
    public ArrayList<ArrayList<Integer>> pathSum(TreeNode root, int sum) {
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        if(root == null)	return result;
        ArrayList<Integer> cur = new ArrayList<Integer>();
        pathSumRec(root,sum,cur,result);
        return result;
    }

    //permutation sequence
    //2^(n-1) repeat
    public String getPermutation(int n, int k) {
    	int []tmp = new int[n+1];
    	tmp[0] = tmp[1] = 1;
    	for(int i=2;i<=n;i++)
    		tmp[i] = i*tmp[i-1];
    	boolean []flag = new boolean[n+1];
    	int []avail = new int[n+1];
    	for(int i=1;i<=n;i++)
    		avail[i] = i;
        String result = "";
        for(int i=n;i>=1;i--){
        	int cur_pos = avail[(k-1)/tmp[i-1]+1];
        	result += cur_pos;
        	flag[cur_pos] = true;
        	int t = 1;
        	k = (k-1) % tmp[i-1] + 1;
        	for(int j=1;j<=n; j++){
        		if(flag[j] == false)
        			avail[t++] = j;
        	}
        	
        }
        return result;
    }
    
    //permutations
    public void permuteRec(int []num,int n,boolean []record, int []cur,ArrayList<ArrayList<Integer>> result){
    	if(n==num.length){
    		ArrayList<Integer> tmp_res = new ArrayList<Integer>();
    		for(int i=0;i<cur.length;i++)
    			tmp_res.add(cur[i]);
    		result.add(tmp_res);
    		return;
    	}
    	for(int i=0;i<num.length;i++){
    		if(record[i] == false){
    			record[i] = true;
    			cur[n] = num[i];
    			permuteRec(num,n+1,record,cur,result);
    			record[i] = false;
    		}
    	}
    }
    public ArrayList<ArrayList<Integer>> permute(int[] num) {
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
    	if(num.length == 0)	return result;
    	Arrays.sort(num);
    	boolean []record = new boolean[num.length];
    	int []cur = new int[num.length];
    	permuteRec(num,0,record,cur,result);
        return result;
    }

    //permutations 2
    public void permuteUniRec(int []num, int n, boolean []record, int []cur, ArrayList<ArrayList<Integer>> result){
    	if(n==num.length){
    		ArrayList<Integer> tmp_res = new ArrayList<Integer>();
    		for(int i=0; i<cur.length; i++)
    			tmp_res.add(cur[i]);
    		result.add(tmp_res);
    		return;
    	}
    	int prev = 0;
    	int i=0;
    	for(i=0; i<num.length; i++){
    		if(record[i] == false){
    			prev = num[i];
    			record[i] = true;
    			cur[n] = num[i];
    			permuteUniRec(num,n+1,record,cur,result);
    			record[i] = false;
    			break;
    		}
    	}
    	for(;i<num.length;i++){
    		if(record[i] == true || num[i] == prev)	continue;
    		prev = num[i];
    		record[i] = true;
    		cur[n] = num[i];
    		permuteUniRec(num,n+1,record,cur,result);
    		record[i] = false;
    	}
    }
    public ArrayList<ArrayList<Integer>> permuteUnique(int[] num) {
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        if(num.length == 0)	return result;
        Arrays.sort(num);
        boolean []flag = new boolean[num.length];
        int []cur = new int[num.length];
        permuteUniRec(num,0,flag,cur,result);
        return result;
    }
    
    //plus one
    public int[] plusOne(int[] digits) {
        int []tmp = new int[digits.length];
        tmp[digits.length-1] = (digits[digits.length-1]+1)%10;
        int carry = (digits[digits.length-1]+1)/10;
        for(int i=digits.length-2; i>=0; i--){
        	tmp[i] = (digits[i]+carry)%10;
        	carry = (digits[i]+carry)/10;
        }
        if(carry == 0)	return tmp;
        int []result = new int[digits.length+1];
        for(int i=1;i<=digits.length;i++)
        	result[i] = tmp[i-1];
        result[0] = carry;
        return result;
    }

    //populating next right pointers in each node
    public void connect(TreeLinkNode root) {
        // Start typing your Java solution below
        // DO NOT write main() function
        if(root == null)	return;
        TreeLinkNode tmp = root.next;
        TreeLinkNode leftmostchild = null;
        while(tmp != null){
        	if(tmp.left != null){
        		leftmostchild = tmp.left;
        		break;
        	}
        	if(tmp.right != null){
        		leftmostchild = tmp.right;
        		break;
        	}
        	tmp = tmp.next;
        }
        if(root.right != null)
        	root.right.next = leftmostchild;
        if(root.left != null){
        	if(root.right != null)
        		root.left.next = root.right;
        	else
        		root.left.next = leftmostchild;
        }
        connect(root.right);
        connect(root.left);
    }

    //pow(x,n)
    public double pow(double x, int n) {
        if(n==0)	return 1;
        if(n==1)	return x;
        if(n==-1)	return 1/x;
        if(n%2 == 0){
        	double tmp = pow(x,n/2);
        	return tmp * tmp;
        }else{
        	double tmp = pow(x,n/2);
        	double tmp2 = n >= 0 ? x : 1/x;
        	return tmp * tmp * tmp2;
        }  
    }
    
    //recover binary tree
    TreeNode prev = null;
    public void findNodes(TreeNode root, ArrayList<TreeNode> nodes){
    	if(root == null)	return ;
    	findNodes(root.left,nodes);
    	if(prev != null && root.val < prev.val){
    		nodes.add(prev);
    		nodes.add(root);
    		if(nodes.size() == 4)
    			return;
    	}
    	prev = root;
    	findNodes(root.right, nodes);
    }
    public void recoverTree(TreeNode root) {
        ArrayList<TreeNode> nodes = new ArrayList<TreeNode>();
        prev = null;
        findNodes(root,nodes);
        TreeNode n1 = null, n2 = null;
        if(nodes.size() == 2){
        	n1 = nodes.get(0);
        	n2 = nodes.get(1);
        }else if(nodes.size() == 4){
        	n1 = nodes.get(0);
        	n2 = nodes.get(3);
        }
        int tmp = n1.val;
        n1.val = n2.val;
        n2.val = tmp;
    }

    //remove duplicates from sorted array
    public int removeDuplicates(int[] A) {
        int result = 0;
        if(A.length == 0)	return result;
        if(A.length == 1)	return 1;
        int prev = A[0], cur_pos = 1;
        for(int i=1; i<A.length; i++){
        	if(A[i] == prev)	continue;
        	A[cur_pos++] = A[i];
        	prev = A[i];
        }
        result = cur_pos;
        return result;
    }
    
    //remove duplicates from sorted array 2
    public int removeDuplicates2(int[] A) {
        if(A.length == 0 || A.length == 1 || A.length == 2)		return A.length;
        int result = 0;
        int prev = A[0], cnt = 1, cur_index = 1;
        for(int i=1;i<A.length;i++){
        	if(A[i] == prev && cnt >= 2)	continue;
        	A[cur_index++] = A[i];
        	if(A[i] == prev)
        		cnt ++;
        	else{
        		prev = A[i];
        		cnt = 1;
        	}
        }
        result = cur_index;
        return result;
    }
    
    //remove duplicates from sorted list
    public ListNode deleteDuplicates(ListNode head) {
        if(head == null)	return null;
        ListNode result = head;
        ListNode cur = result;
        int prev = head.val;
        cur = head;
        head = head.next;
        cur.next = null;
        while(head != null){
        	if(head.val == prev){
        		ListNode tmp = head;
        		head = head.next;
        		tmp.next = null;
        		continue;
        	}
        	prev = head.val;
        	cur.next = head;
        	cur = head;
        	head = head.next;
        	cur.next = null;
        }
        return result;
    }

    //remove duplicates from sorted list 2
    public ListNode deleteDuplicates2(ListNode head) {
        ListNode result = null;
        if(head == null)	return result;
        ListNode cur = head.next;
        ListNode prev = head;
        ListNode rHead = null;
        boolean isRep = false;
        while(cur != null){
        	if(cur.val == prev.val)
        		isRep = true;
        	else{
        		if(!isRep){
        			if(result == null){
        				result = rHead = prev;
        				result.next = null;
        			}else{
        				result.next = prev;
        				result = prev;
        				result.next = null;
        			}
        		}
        		isRep = false;
        		prev = cur;
        	}
        	cur = cur.next;
        }
        if(!isRep){
        	if(result == null)
        		result = rHead = prev;
        	else{
        		result.next = prev;
        		result = prev;
        	}
        	result.next = null;
        }
        
        return rHead;
    }

    //remove element
    public int removeElement(int[] A, int elem) {
        int cur_pos = 0;
        for(int i=0;i<A.length; i++){
        	if(A[i] != elem)
        		A[cur_pos++] = A[i];
        }
        return cur_pos;
    }

    //remove Nth node from end of list
    public ListNode removeNthFromEnd(ListNode head, int n) {
        ListNode result = null;
        if(head == null)	return result;
        ListNode first = head, second = head, prev = null;
        for(int i=1; i<=n; i++)		second = second.next;
        while(second != null){
        	prev = first;
        	first = first.next;
        	second = second.next;
        }
        if(prev == null){
        	result = first.next;
        }else{
        	result = head;
        	prev.next = first.next;
        }
        return result;
    }

    //restore IP address
    public void findValidIPRec(String s, int n, int start, char[]record, ArrayList<String> result){
    	if(n == 4 && start == s.length()){
    		String tmp_s = new String(record);
    		result.add(tmp_s);
    		return;
    	}
    	if(n == 4 || start >= s.length())	return;
    	if(n == 3){
    		if(s.charAt(start) == '0'){
    			if(start == s.length()-1){
    				record[start+3] = '0';
    				findValidIPRec(s,n+1,start+1,record,result);
    			}
    		}else{
        		int tmp = 0, cnt = 0;
        		for(int i=start;i<s.length() && cnt < 3;i++, cnt++)
        			tmp = tmp * 10 + s.charAt(i) - '0';
        		if(tmp > 255)	return;
        		for(int i=start+3;i<start+3+cnt;i++)
        			record[i] = s.charAt(i-3);
        		findValidIPRec(s,n+1,start+cnt,record,result);
    		}
    	}else{
    		if(s.charAt(start) == '0'){
    			record[start+n] = '0';
    			record[start+n+1] = '.';
    			findValidIPRec(s,n+1,start+1,record,result);
    		}else{
    			char c1 = s.charAt(start);
    			record[start+n] = c1; record[start+n+1] = '.';
    			findValidIPRec(s,n+1,start+1,record,result);
    			if(start+1 <s.length()){
    				char c2 = s.charAt(start+1);
    				record[start+n] = c1; record[start+n+1] = c2; record[start+n+2] = '.';
    				findValidIPRec(s,n+1,start+2,record,result);
    				if(start+2 < s.length()){
    					char c3 = s.charAt(start+2);
    					int tmp = (c1-'0')*100+(c2-'0')*10+c3-'0';
    					if(tmp <= 255){
    						record[start+n] = c1; record[start+n+1] = c2; record[start+n+2] = c3; record[start+n+3] = '.';
    						findValidIPRec(s,n+1,start+3,record,result);
    					}
    				}
    			}
    		}
    	}
    }
    public ArrayList<String> restoreIpAddresses(String s) {
        ArrayList<String> result = new ArrayList<String>();
        char []record = new char[s.length()+3];
        findValidIPRec(s,0,0,record,result);
        return result;
    }

    //reverse integer
    public int reverse(int x) {
        boolean isNegative = x < 0 ? true : false;
        if(isNegative) x = -x;	//may overflow
        int tmp = x, high_base = 1, result = 0, pos_num = 1;
        while(tmp != 0){
        	tmp /= 10;
        	high_base *= 10;
        }
        high_base /= 10;
        while(x != 0){
        	int ori_low = x%10;
        	int ori_high = x/high_base;
        	if(high_base == 1)
        		result += ori_low * pos_num;
        	else
        		result += (ori_low * high_base + ori_high)*pos_num;
        	x %= high_base;
        	x /= 10;
        	high_base /= 100;
        	pos_num *= 10;
        }
        result = isNegative ? -result : result;
        return result;
    }

    //reverse linked list 2
    public ListNode reverseBetween(ListNode head, int m, int n) {
        if(head == null)	return null;
        ListNode cur = head, prev = null, prev_prev = null, start = null;
        int cnt = 1;
        while(cnt <= m){
        	prev_prev = prev;
        	prev = cur;
        	cur = cur.next;
        	cnt ++;
        }
        start = prev;
        while(cnt <= n){
        	ListNode tmp = cur.next;
        	cur.next = prev;
        	prev = cur;
        	cur = tmp;
        	cnt ++;
        }
        if(prev_prev != null)
        	prev_prev.next = prev;
        start.next = cur;
        if(m == 1)	return prev;
        return head;
    }

    //reverse nodes in k-group
    //consider a simple and straight-forward solution
    public ListNode reverseKGroup(ListNode head, int k) {
        if(head == null)	return null;
        if(k==1)	return head;
        int num = 0, group_num = 0;
        ListNode cur = head, prev = null,group_prev1 = null, group_prev2 = null;
        while(cur != null){
        	cur = cur.next;
        	num ++;
        }
        cur = head;
        while(cur != null){
        	if(group_num == 0){
        		if(group_prev1 != null)
        			group_prev2 = group_prev1;
        		if(num >= k){
        			num --;
        			group_prev1 = cur;
        			group_num = (group_num+1)%k;
        			ListNode tmp = cur.next;
        			prev = cur;
        			cur.next = null;
        			cur = tmp;
        			continue;
        		}else{
        			break;
        		}
        	}
        	if(group_num == k-1){
        		if(group_prev2 != null)
        			group_prev2.next = cur;
        		else
        			head = cur;
        	}
        	ListNode tmp = cur.next;
        	cur.next = prev;
        	prev = cur;
        	cur = tmp;
        	num --;
        	group_num = (group_num+1)%k;
        }
        if(cur != null && group_prev1 != null)
        	group_prev1.next = cur;
        return head;
    }

    //rotate image
    public void rotate(int[][] matrix) {
        for(int d=0; d<matrix.length/2; d++){
        	for(int i=0;d+i<matrix.length-d-1;i++){
        		int tmp = matrix[d][d+i];	
        		matrix[d][d+i] = matrix[matrix.length-1-d-i][d];
        		matrix[matrix.length-1-d-i][d] = matrix[matrix.length-1-d][matrix.length-1-d-i];
        		matrix[matrix.length-1-d][matrix.length-1-d-i] = matrix[d+i][matrix.length-1-d];
        		matrix[d+i][matrix.length-1-d] = tmp;
        	}
        }
    }
    
    //rotate list
    public ListNode rotateRight(ListNode head, int n) {
        ListNode cur = head, fast = head, prev = null;
        if(head == null || n == 0)	return head;
        int num = 0;
        while(cur != null){
        	cur = cur.next;
        	num ++;
        }
        cur = head;
        n = n % num;
        if(n == 0)	return head;
        for(int i=1;i<n && fast != null;i++)
        	fast = fast.next;
        if(fast == null)	return head;
        while(fast.next != null){
        	prev = cur;
        	cur = cur.next;
        	fast = fast.next;
        }
        if(prev != null){
        	prev.next = null;
        	fast.next = head;
        	head = cur;
        }
        return head;
    }

    //same tree
    public boolean isSameTree(TreeNode p, TreeNode q) {
        if(p == null && q == null)	return true;
        if(p == null || q == null)	return false;
        boolean left = isSameTree(p.left,q.left);
        boolean right = isSameTree(p.right,q.right);
        return left && right && (p.val == q.val);
    }
    
    //search a 2d matrix
    public boolean searchMatrix(int[][] matrix, int target) {
        int sr = 0, er = matrix.length-1, tr = -1;
        while(sr <= er){
        	int mid = (sr+er)/2;
        	if(target >= matrix[mid][0] && target <= matrix[mid][matrix[0].length-1]){
        		tr = mid;
        		break;
        	}
        	else if(target < matrix[mid][0])
        		er = mid-1;
        	else
        		sr = mid+1;
        }
        if(tr == -1)	return false;
        sr = 0; er = matrix[0].length-1;
        while(sr <= er){
        	int mid = (sr+er)/2;
        	if(matrix[tr][mid] == target)	return true;
        	if(matrix[tr][mid] < target)	sr = mid + 1;
        	else	er = mid - 1;
        }
        return false;
    }

    //search for a range
    public int[] searchRange(int[] A, int target) {
        int s = 0, e = A.length-1;
        int []result = {-1,-1};
        while(s <= e){
        	int mid = (s+e)/2;
        	if(target < A[mid]) e--;
        	else if(target > A[mid]) s++;
        	else{
        		if(mid == 0 || A[mid-1] < target){
        			result[0] = mid;
        			break;
        		}else{
        			e --;
        		}
        	}
        }
        s = 0; e = A.length-1;
        while(s <= e){
        	int mid = (s+e)/2;
        	if(target < A[mid]) e--;
        	else if(target > A[mid]) s++;
        	else{
        		if(mid == A.length-1 || A[mid+1] > target){
        			result[1] = mid;
        			break;
        		}else{
        			s++;
        		}
        	}
        }
        return result;
    }

    //search in rotated array list
    public int searchRec(int []A, int target, int s, int e){
    	if(s > e)	return -1;
    	int mid = (s+e)/2;
    	if(A[mid] == target)	return mid;
    	if(A[mid] < A[s]){
    		if(target > A[mid]){
    			int tmp_res = searchRec(A,target,s,mid-1);
    			if(tmp_res != -1)
    				return tmp_res;
    			else
    				return searchRec(A,target,mid+1,e);
    		}
    		else
    			return searchRec(A,target,s,mid-1);
    	}else if(A[mid] > A[e]){
    		if(target > A[mid])
    			return searchRec(A,target,mid+1,e);
    		else{
    			int tmp_res = searchRec(A,target,mid+1,e);
    			if(tmp_res != -1)
    				return tmp_res;
    			else
    				return searchRec(A,target,s,mid-1);
    		}
    	}else{
    		return target > A[mid] ? searchRec(A,target,mid+1,e) : searchRec(A,target,s,mid-1);
    	}
    }
    public int search(int[] A, int target) {
        return searchRec(A,target,0,A.length-1);
    }

    //search in rotated array list 2
    public boolean searchRec2(int []A, int target, int s, int e){
    	if(s>e)		return false;
    	int mid = (s+e)/2;
    	if(A[mid] == target)	return true;
    	if(A[mid] > target){
    		if(A[mid] > A[e]){
    			if(searchRec2(A,target,s,mid-1) == true)
    				return true;
    			return searchRec2(A,target,mid+1,e);
    		}else if(A[mid] < A[s]){
    			return searchRec2(A,target,s,mid-1);
    		}else if(A[mid] > A[s] && A[mid] < A[e]){
    			return searchRec2(A,target,s,mid-1);
    		}else{
    			if(searchRec2(A,target,s,mid-1) == true)
    				return true;
    			return searchRec2(A,target,mid+1,e);
    		}
    	}else{
    		if(A[mid] > A[e]){
    			return searchRec2(A,target,mid+1,e);
    		}else if(A[mid] < A[s]){
    			if(searchRec2(A,target,mid+1,e) == true)
    				return true;
    			return searchRec2(A,target,s,mid-1);
    		}else if(A[mid] > A[s] && A[mid] < A[e]){
    			return searchRec2(A,target,mid+1,e);
    		}else{
    			if(searchRec2(A,target,s,mid-1) == true)
    				return true;
    			return searchRec2(A,target,mid+1,e);
    		}
    	}
    }
    public boolean search2(int[] A, int target) {
        return searchRec2(A,target,0,A.length-1);
    }
    
    //search insert position
    public int searchInsert(int[] A, int target) {
        if(A.length == 0)	return 0;
        for(int i=0;i<A.length;i++){
        	if(A[i] >= target)	return i;
        }
        return A.length;
    }

    //set matrix zeros
    public void setZeroes(int[][] matrix) {
        int r = -1 ,c = -1;
        if(matrix.length == 0 || matrix[0].length == 0)		return;
        for(int i=0;i<matrix.length;i++){
        	for(int j=0;j<matrix[0].length;j++){
        		if(matrix[i][j] == 0){
        			if(c == -1 && r == -1){
        				r = i;
        				c = j;
        			}else{
        				matrix[r][j] = 0;
        				matrix[i][c] = 0;
        			}
        		}
        	}
        }
        if(r == -1 && c == -1) 		return;
        for(int i=0;i<matrix[0].length;i++){
        	if(matrix[r][i] == 0){
        		if(i == c)	continue;
        		for(int j=0;j<matrix.length;j++)
        			matrix[j][i] = 0;
        	}
        }
        for(int i=0; i<matrix.length; i++){
        	if(matrix[i][c] == 0){
        		for(int j=0;j<matrix[0].length;j++)
        			matrix[i][j] = 0;
        	}
        }
        for(int i=0; i<matrix.length; i++)
        	matrix[i][c] = 0;
    }

    //simplify path
    public String simplifyPath(String path) {
        if(path == null || path.length() == 0)	return "";
        String result = "";
        String []strs = path.split("/");
        Stack<String> st = new Stack<String>();
        for(int i=0;i<strs.length;i++){
        	if(strs[i].equals("") || strs[i].equals(".")) continue;
        	if(strs[i].equals("..")){
        		if(!st.isEmpty())
        			st.pop();
        	}else{
        		st.push(strs[i]);
        	}
        }
        if(st.isEmpty())
        	return "/";
        String []tmp = new String[strs.length];
        int cur = 0;
        while(!st.isEmpty()){
        	tmp[cur++] = st.pop();
        }
        cur --;
        for(;cur >= 0; cur--)
        	result += "/" + tmp[cur];
        return result;
    }

    //sort colors
    public void sortColors(int[] A) {
        if(A.length == 0)	return;
        int s = 0, e = A.length-1;
        for(int i=0;i<=e; i++){
        	if(A[i] == 2){
        		while(A[e] == 2 && e > i)	e--;
        		int tmp = A[e];
        		A[e] = A[i];
        		A[i] = tmp;
        		e --;
        		if(A[i] == 0){
        			tmp = A[s];
        			A[s] = A[i];
        			A[i] = tmp;
        			s ++;
        		}
        	}else if(A[i] == 0){
        		int tmp = A[s];
        		A[s] = A[i];
        		A[i] = tmp;
        		s ++;
        	}
        }
    }
    
    //spiral matrix
    public ArrayList<Integer> spiralOrder(int[][] matrix) {
        ArrayList<Integer> result = new ArrayList<Integer> ();
        if(matrix.length == 0 || matrix[0].length == 0)		return result;
        int short_size = matrix.length < matrix[0].length ? matrix.length : matrix[0].length;
        int d = 0;
        for(d=0;d<short_size/2; d++){
        	for(int i=d;i<matrix[0].length-d;i++)
        		result.add(matrix[d][i]);
        	for(int i=d+1; i<matrix.length-d; i++)
        		result.add(matrix[i][matrix[0].length-d-1]);
        	for(int i=matrix[0].length-d-2;i>=d;i--)
        		result.add(matrix[matrix.length-d-1][i]);
        	for(int i=matrix.length-d-2;i>d;i--)
        		result.add(matrix[i][d]);
        }
        if(short_size%2 == 1){
        	if(short_size == matrix.length){
        		for(int i=d; i<matrix[0].length-d;i++)
        			result.add(matrix[d][i]);
        	}else{
        		for(int i=d;i<matrix.length-d;i++)
        			result.add(matrix[i][d]);
        	}
        }
        return result;
    }
    
    //spiral matrix 2
    public int[][] generateMatrix(int n) {
        int [][]result = new int[n][n];
        if(n<=0)	return result;
        int cur = 1;
        for(int d=0; d < n/2; d++){
        	for(int i=d;i<n-d;i++)
        		result[d][i] = cur ++;
        	for(int i=d+1;i<n-d;i++)
        		result[i][n-d-1] = cur ++;
        	for(int i=n-d-2;i>=d;i--)
        		result[n-d-1][i] = cur ++;
        	for(int i=n-d-2;i>d;i--)
        		result[i][d] = cur ++;
        }
        if(n%2 == 1)
        	result[n/2][n/2] = n*n;
        return result;
    }

    //sqrt(x)
    //long should avoid overflow
    //try newton mwthod someday
    public int sqrt(int x) {
    	if(x <= 0)	return 0;
    	if(x == 1)	return 1;
    	long lx = x;
    	long s = 1, e = x;
    	while(s<=e){
    		long mid = (s+e)/2;
    		long tmp = mid * mid;
    		if(tmp == lx)	return (int)mid;
    		if(tmp > x)		e = mid - 1;
    		else 	s = mid + 1;
    	}
    	return (int)e;
    }

    //string to integer
    //+,-,' ' before the num
    public int atoi(String str) {
        int result = 0;
        if(str == null || str.length() == 0)	return result;
        boolean isFound = false, isNeg = false;
        long tmp_result = 0;
        for(int i=0;i<str.length();i++){
        	if(str.charAt(i) != '-' && str.charAt(i) != '+' && !(str.charAt(i) >= '0' && str.charAt(i) <= '9')){
        		if(!isFound && str.charAt(i) == ' ')	continue;
        		if(!isFound)	return 0;
        		return isNeg ? (int)-tmp_result : (int)tmp_result;
        	}
        	if(str.charAt(i) == '-'){
        		if(isFound)
        			return isNeg ? (int)-tmp_result : (int)tmp_result;
        		isNeg = isFound = true;
        		continue;
        	}
        	if(str.charAt(i) == '+'){
        		if(isFound)
        			return isNeg ? (int)-tmp_result : (int)tmp_result;
        		isNeg = false;
        		isFound = true;
        		continue;
        	}
        	isFound = true;
        	tmp_result = tmp_result * 10 + (str.charAt(i) - '0');
        	if(!isNeg && tmp_result >= Integer.MAX_VALUE)
        		return Integer.MAX_VALUE;
        	if(isNeg && (-tmp_result) <= Integer.MIN_VALUE)
        		return Integer.MIN_VALUE;
        }
        result = isNeg ? (int)-tmp_result : (int)tmp_result;
        return result;
    }
    
    //subsets
    public void subsetRec(int []S, int n, int num, ArrayList<ArrayList<Integer>> result, int[]record){
    	if(n == S.length){
    		ArrayList<Integer> tmp_res = new ArrayList<Integer>();
    		for(int i=0;i<num;i++)
    			tmp_res.add(record[i]);
    		result.add(tmp_res);
    		return;
    	}
    	subsetRec(S,n+1,num,result,record);
    	record[num] = S[n];
    	subsetRec(S,n+1,num+1,result,record);
    }
    public ArrayList<ArrayList<Integer>> subsets(int[] S) {
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        if(S.length == 0){
            ArrayList<Integer> empty = new ArrayList<Integer>();
            result.add(empty);
        	return result;
        }
        Arrays.sort(S);
        int []record = new int[S.length];
        subsetRec(S,0,0,result,record);
        return result;
    }

    //subsets 2
    //choosing empty : skip the same ones
    public void subsetDupRec(int []S, int n, int num, ArrayList<ArrayList<Integer>> result, int []record){
    	if(n == S.length){
    		ArrayList<Integer> tmp_res = new ArrayList<Integer>();
    		for(int i=0;i<num; i++)
    			tmp_res.add(record[i]);
    		result.add(tmp_res);
    		return;
    	}
    	record[num] = S[n];
    	subsetDupRec(S,n+1,num+1,result,record);
    	int prev = S[n];
    	for(;n<S.length && S[n] == prev;n++);
    	subsetDupRec(S,n,num,result,record);
    }
    public ArrayList<ArrayList<Integer>> subsetsWithDup(int[] num) {
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        if(num.length == 0){
        	result.add(new ArrayList<Integer>());
        	return result;
        }
        Arrays.sort(num);
        int []record = new int[num.length];
        subsetDupRec(num,0,0,result,record);
        return result;
    }
    
    //substring with concatenation of all words
    public ArrayList<Integer> findSubstring(String S, String[] L) {
        ArrayList<Integer> result = new ArrayList<Integer>();
        if(S.length() == 0 || L.length == 0 || L[0].length() == 0)	return result;
        int word_len = L[0].length(), prev_start = 0, left = L.length;
        HashMap<String,Integer> dict = new HashMap<String,Integer>();
        HashMap<String,Integer> hm = new HashMap<String,Integer>();
        for(int i=0;i<L.length; i++){
        	if(!dict.containsKey(L[i]))
        		dict.put(L[i], 1);
        	else
        		dict.put(L[i], dict.get(L[i])+1);
        }
        int i = 0;
        for(i=0;i<S.length();i++){
        	if(i+L.length*L[0].length()-1 >= S.length())	break;
        	hm.clear();
        	left = L.length;
        	for(int j=i;j<i+L.length*L[0].length();j+=word_len){
        		String cur_word = S.substring(j,j+word_len);
        		if(!dict.containsKey(cur_word) || (hm.containsKey(cur_word) && hm.get(cur_word) >= dict.get(cur_word)))
        			break;
        		if(!hm.containsKey(cur_word))
        			hm.put(cur_word, 1);
        		else
        			hm.put(cur_word, hm.get(cur_word)+1);
        		left --;
        	}
        	if(left == 0)
        		result.add(i);
        }
        return result;
    }

    //swap nodes in pairs
    public ListNode swapPairs(ListNode head) {
        ListNode result = head;
        if(head == null || head.next == null)	return result;
        ListNode cur = head, prev = null, prevprev = null;
        result = cur.next;
        int num = 0;
        while(cur != null){
        	if(num == 0){
        		if(cur.next == null && prevprev != null)
        			prevprev.next = cur;
        		prev = cur;
        		ListNode t_next = cur.next;
        		cur.next = null;
        		cur = t_next;
        		num = (num + 1) % 2;
        	}else{
        		if(prevprev != null)
        			prevprev.next = cur;
        		prevprev = prev;
        		ListNode t_next = cur.next;
        		cur.next = prev;
        		prev = cur;
        		cur = t_next;
        		num = (num+1)%2;
        	}
        }
        return result;
    }

    //symmetric tree
    public boolean isSymmetric(TreeNode root1, TreeNode root2){
    	if(root1 == null && root2 == null)	return true;
    	if(root1 == null || root2 == null)	return false;
    	if(root1.val != root2.val) return false;
    	return isSymmetric(root1.right,root2.left) && isSymmetric(root1.left,root2.right);
    }
    public boolean isSymmetric(TreeNode root) {
    	if(root == null)	return true;
        return isSymmetric(root.left, root.right);
    }

    //text justification
    public ArrayList<String> fullJustify(String[] words, int L) {
        ArrayList<String> result = new ArrayList<String>();
        if(words.length == 0)	return result;
        int cur_len = 0;
        ArrayList<String> tmp = new ArrayList<String>();
        for(int i=0; i<words.length; i++){
        	String cur_word = words[i];
        	if(cur_len == 0){
        		cur_len = cur_word.length();
        		tmp.add(cur_word);
        		continue;
        	}
        	if(cur_len + cur_word.length() + 1 > L){
        		int space_num = L - cur_len;
        		int should_space_num = tmp.size()-1;
        		if(should_space_num == 0){
        			String line = tmp.get(0);
        			for(int j=line.length()+1;j<=L;j++)
        				line = line + " ";
        			result.add(line);
        		}else{
            		int space_per_word = space_num/should_space_num;
            		int word_num_with_add_space = space_num % should_space_num;
            		String line = "";
            		for(int j=0;j<tmp.size()-1;j++){
            			line = line + tmp.get(j);
            			for(int k=1;k<=space_per_word; k++)
            				line = line + " ";
            			if(word_num_with_add_space > 0)
            				line = line + " ";
            			word_num_with_add_space --;
            		}
            		line = line + tmp.get(tmp.size()-1);
            		result.add(line);	
        		}
        		tmp.clear();
        		tmp.add(cur_word);
        		cur_len = cur_word.length();
        	}else{
        		cur_len += cur_word.length()+1;
        		cur_word = " " + cur_word;
        		tmp.add(cur_word);
        	}
        }
        if(tmp.size() != 0){
        	String line = "";
        	for(int i=0;i<tmp.size();i++){
        		line = line + tmp.get(i);
        	}
        	for(int i=line.length()+1;i<=L;i++)
        		line = line + " ";
        	result.add(line);
        }
        return result;
    }

    //trapping rain water
    public int trap(int[] A) {
    	int result = 0;
    	Stack<Integer> x = new Stack<Integer>();
    	Stack<Integer> h = new Stack<Integer>();
    	for(int i=0;i<A.length;i++){
    		if(x.isEmpty()){
    			x.push(i);
    			h.push(A[i]);
    		}else{
    			if(A[i] <= h.peek()){
    				x.push(i);
    				h.push(A[i]);
    			}else{
    				int bottom = h.pop();
    				x.pop();
    				while(!x.isEmpty() && h.peek() <= A[i]){
    					result += (h.peek()-bottom)*(i-1-x.pop());
    					bottom = h.pop();
    				}
    				if(!x.isEmpty()){
    					result += (A[i]-bottom)*(i-1-x.peek());
    				}
    				x.push(i);
    				h.push(A[i]);
    			}
    		}
    	}
        return result;
    }

    //triangle
    public int minimumTotal(ArrayList<ArrayList<Integer>> triangle) {
        int result = 0;
        if(triangle == null || triangle.size() == 0)	return result;
        int []prev = new int[triangle.size()];
        int []cur = new int[triangle.size()];
        prev[0] = triangle.get(0).get(0);
        for(int i=1;i<triangle.size();i++){
        	for(int j=0;j<=i;j++){
        		if(j == 0)	cur[j] = prev[0] + triangle.get(i).get(j);
        		else if(j == i)	cur[j] = prev[i-1] + triangle.get(i).get(j);
        		else cur[j] = prev[j-1] < prev[j] ? prev[j-1] + triangle.get(i).get(j) : prev[j] + triangle.get(i).get(j);
        	}
        	for(int j=0;j<=i;j++)
        		prev[j] = cur[j];
        }
        result = prev[0];
        for(int i=0;i<triangle.size();i++)
        	result = prev[i] < result ? prev[i] : result;
        return result;
    }

    //two sum
    public int[] twoSum(int[] numbers, int target) {
        int []result = new int[2];
        HashMap<Integer,Integer> hm = new HashMap<Integer,Integer>();
        for(int i=0;i<numbers.length;i++){
        	hm.put(numbers[i], i);
        }
        for(int i=0;i<numbers.length;i++){
        	int rest = target-numbers[i];
        	if(hm.containsKey(rest)){
        		result[0] = i+1;
        		result[1] = hm.get(rest)+1;
        		break;
        	}
        }
        return result;
    }

    //unique binary trees
    public int numTrees(int n) {
        if(n == 0)	return 0;
        if(n == 1)	return 1;
        if(n == 2)	return 2;
        int []tmp = new int[n+1];
        tmp[0] = 0; tmp[1] = 1; tmp[2] = 2;
        for(int i=3;i<=n;i++){
        	tmp[i] = 2*tmp[i-1];
        	for(int j=1;j<i-1;j++)
        		tmp[i] += tmp[j] * tmp[i-1-j];
        }
        return tmp[n];
    }

    //unique binary trees 2
    ArrayList<TreeNode> generateTreeRec(int s, int e){
    	ArrayList<TreeNode> result = new ArrayList<TreeNode>();
    	if(s > e){
    		result.add(null);
    		return result;
    	}
    	if(s == e){
    		TreeNode n = new TreeNode(s);
    		result.add(n);
    		return result;
    	}
    	for(int i=s;i<=e;i++){
    		ArrayList<TreeNode> leftT = generateTreeRec(s,i-1);
    		ArrayList<TreeNode> rightT = generateTreeRec(i+1,e);
    		for(int j1 = 0; j1<leftT.size(); j1++){
    			for(int j2 = 0; j2<rightT.size(); j2++){
    				TreeNode cur = new TreeNode(i);
    				cur.left = leftT.get(j1);
    				cur.right = rightT.get(j2);
    				result.add(cur);
    			}
    		}
    	}
    	return result;
    }
    public ArrayList<TreeNode> generateTrees(int n) {
        return generateTreeRec(1,n);
    }

    //unique paths
    public int uniquePaths(int m, int n) {
    	if(m <=0 || n <= 0)		return 0;
    	int [][]record = new int[m][n];
    	for(int i=0;i<m;i++)	record[i][0] = 1;
    	for(int i=0;i<n;i++)	record[0][i] = 1;
        for(int i=1;i<m;i++)
        	for(int j=1;j<n;j++)
        		record[i][j] = record[i][j-1] + record[i-1][j];
        return record[m-1][n-1];
    }

    //unique paths 2
    public int uniquePathsWithObstacles(int[][] obstacleGrid) {
        if(obstacleGrid.length == 0 || obstacleGrid[0].length == 0)
        	return 0;
        int [][]record = new int[obstacleGrid.length][obstacleGrid[0].length];
        boolean flag = false;
        for(int i=0;i<obstacleGrid.length;i++){
        	if(obstacleGrid[i][0] == 1)		flag = true;
        	record[i][0] = flag ? 0 : 1;
        }
        flag = false;
        for(int i=0;i<obstacleGrid[0].length;i++){
        	if(obstacleGrid[0][i] == 1)		flag = true;
        	record[0][i] = flag ? 0 : 1;
        }
        for(int i=1;i<obstacleGrid.length;i++)
        	for(int j=1;j<obstacleGrid[0].length;j++){
        		if(obstacleGrid[i][j] == 1)
        			record[i][j] = 0;
        		else
        			record[i][j] = record[i-1][j] + record[i][j-1];
        	}
        return record[obstacleGrid.length-1][obstacleGrid[0].length-1];
    }

    //valid number
    public boolean isNum(String s){
    	for(int i=0;i<s.length();i++){
    		if(!(s.charAt(i) >= '0' && s.charAt(i) <= '9'))
    			return false;
    	}
    	return true;
    }
    public boolean isNum2(String s){
    	if(s.length() == 0)		return false;
    	if(s.contains(".")){
    		int pos = s.indexOf(".");
    		if(s.length() == 1)		return false;
    		return isNum(s.substring(0,pos)) && isNum(s.substring(pos+1,s.length()));
    	}else
    		return isNum(s);
    }
    public boolean isNumber(String s) {
    	int backward_space = -1;
    	for(int i=s.length()-1;i>=0;i--){
    		if(s.charAt(i) != ' ')	{
    			backward_space = i;
    			break;
    		}
    	}
    	if(backward_space == -1)	return false;
    	s = s.substring(0,backward_space+1);
    	backward_space = s.length();
    	for(int i=0;i<s.length();i++){
    		if(s.charAt(i) != ' '){
    			backward_space = i;
    			break;
    		}
    	}
    	if(backward_space == s.length())	return false;
    	s = s.substring(backward_space);
    	if(s.contains("e")){
    		if(s.charAt(0) == '+' || s.charAt(0) == '-')
    			s = s.substring(1,s.length());
    		int pos = s.indexOf("e");
    		if(pos == s.length()-1)		return false;
    		if(s.charAt(pos+1) == '+' || s.charAt(pos+1) == '-'){
    			if(pos+2 == s.length())		return false;
    			return isNum2(s.substring(0,pos)) && isNum(s.substring(pos+2,s.length()));
    		}
    		else
    			return isNum2(s.substring(0,pos)) && isNum(s.substring(pos+1,s.length()));		
    	}else{
    		if(s.charAt(0) == '+' || s.charAt(0) == '-')
    			return isNum2(s.substring(1,s.length()));
    		else
    			return isNum2(s);
    	}
    }

    //valid palindrome
    public boolean isPalindrome(String s) {
        if(s == null || s.length() == 0)	return true;
        char[]record = new char[s.length()];
        int cur = 0;
        for(int i=0;i<s.length();i++){
        	if((s.charAt(i) >= 'a' && s.charAt(i) <= 'z') || (s.charAt(i) >= '0' && s.charAt(i) <= '9'))
        		record[cur++] = s.charAt(i);
        	else if(s.charAt(i) >= 'A' && s.charAt(i) <= 'Z')
        		record[cur++] = (char)(s.charAt(i) - 'A' + 'a');
        }
        if(cur%2 == 0){
        	for(int i=0;cur/2-1-i>=0 && cur/2+i<cur;i++){
        		if(record[cur/2-1-i] != record[cur/2+i])
        			return false;
        	}
        }else{
        	for(int i=1;cur/2-i >= 0 && cur/2+i < cur;i++){
        		if(record[cur/2-i] != record[cur/2+i])
        			return false;
        	}
        }
        return true;
    }

    //valid parentheses
    public boolean isValid(String s) {
        if(s == null || s.length() == 0)	return true;
        Stack<Character> st = new Stack<Character>();
        for(int i=0;i<s.length();i++){
        	if(s.charAt(i) == '(' || s.charAt(i) == '[' || s.charAt(i) == '{')
        		st.push(s.charAt(i));
        	else{
        		if(st.isEmpty())  return false;
        		char c = st.pop();
        		if(s.charAt(i) == ')'){
        			if(c != '(')
        				return false;
        		}
        		else if(s.charAt(i) == ']'){
        			if(c != '[')
        				return false;
        		}
        		else if(s.charAt(i) == '}'){
        			if(c != '{')
        				return false;
        		}
        	}
        }
        if(st.isEmpty())	return true;
        return false;
    }

    //validate binary search tree
    public int[] isValidBSTRec(TreeNode root){
    	if(root.left == null && root.right == null){
    		int []result = {root.val,root.val};
    		return result;
    	}else if(root.left == null){
    		int []right = isValidBSTRec(root.right);
    		if(right == null)	return null;
    		if(root.val < right[0]){
    			int []result = {root.val,right[1]};
    			return result;
    		}else
    			return null;
    	}else if(root.right == null){
    		int []left = isValidBSTRec(root.left);
    		if(left == null)	return null;
    		if(root.val > left[1]){
    			int []result = {left[0],root.val};
    			return result;
    		}else
    			return null;
    	}else{
    		int []left = isValidBSTRec(root.left);
    		int []right = isValidBSTRec(root.right);
    		if(left == null || right == null)	return null;
    		if(left[1] < root.val && root.val < right[0]){
    			int []result = {left[0],right[1]};
    			return result;
    		}else
    			return null;
    	}
    }
    public boolean isValidBST(TreeNode root) {
    	if(root == null || (root.left == null && root.right == null))	return true;
    	if(root.left == null){
    		int []right = isValidBSTRec(root.right);
    		if(right == null || root.val >= right[0])	return false;
    		return true;
    	}else if(root.right == null){
    		int []left = isValidBSTRec(root.left);
    		if(left == null || root.val <= left[1])		return false;
    		return true;
    	}else{
    		int []left = isValidBSTRec(root.left);
    		int []right = isValidBSTRec(root.right);
    		if(left == null || right == null || left[1] >= root.val || right[0] <= root.val)	return false;
    		return true;
    	}
    }

    //
    int sum = 0;
    public void sumNumRec(TreeNode root, ArrayList<Integer> record){
    	if(root.left == null && root.right == null){
    		int tmp_sum = 0;
    		for(int i=0;i<record.size();i++)
    			tmp_sum = tmp_sum * 10 + record.get(i);
    		tmp_sum = tmp_sum * 10 + root.val;
    		sum += tmp_sum;
    		return;
    	}
    	if(root.left == null){
    		record.add(root.val);
    		sumNumRec(root.right,record);
    		record.remove(record.size()-1);
    	}else if(root.right == null){
    		record.add(root.val);
    		sumNumRec(root.left,record);
    		record.remove(record.size()-1);
    	}else{
    		record.add(root.val);
    		sumNumRec(root.left,record);
    		sumNumRec(root.right,record);
    		record.remove(record.size()-1);
    	}
    }
    public int sumNumbers(TreeNode root) {
    	sum = 0;
        ArrayList<Integer> record = new ArrayList<Integer>();
        if(root == null)	return 0;
        sumNumRec(root,record);
        return sum;
    }

    //Surrounded Regions
    //use Queue(LinkedList) will TLE. Use ArrayDeque instead
    ArrayDeque<int []> q = new ArrayDeque<int []>();
    public boolean isValidPos(char [][]board, boolean [][]record, int x, int y){
    	if(x < 0 || x >= board.length || y < 0 || y >= board[0].length || board[x][y] == 'X' || record[x][y])
    		return false;
    	return true;
    }
    public void solve2(char [][]board, boolean [][]record, int x, int y){
    	int []cur = {x,y};
    	q.offer(cur);
    	record[x][y] = true;
    	while(!q.isEmpty()){
    		int []pos = q.poll();
    		for(int i=-1;i<=1;i+=2){
    			if(isValidPos(board,record,pos[0]+i,pos[1])){
    				int []next = {pos[0]+i,pos[1]};
    				q.offer(next);
    				record[pos[0]+i][pos[1]] = true;
    			}
    			if(isValidPos(board,record,pos[0],pos[1]+i)){
    				int []next = {pos[0],pos[1]+i};
    				q.offer(next);
    				record[pos[0]][pos[1]+i] = true;
    			}
    		}
    	}
    }
    public void solve(char[][] board) {
    	if(board.length == 0 || board[0].length == 0)	return;
        boolean [][]record = new boolean[board.length][board[0].length];
        for(int i=0;i<board.length;i++){
        	if(board[i][0] == 'O' && record[i][0] == false)
        		solve2(board,record,i,0);
        	if(board[i][board[0].length-1] == 'O' && record[i][board[0].length-1] == false)
        		solve2(board,record,i,board[0].length-1);
        }
        for(int i=0;i<board[0].length;i++){
        	if(board[0][i] == 'O' && record[0][i] == false)
        		solve2(board,record,0,i);
        	if(board[board.length-1][i] == 'O' && record[board.length-1][i] == false)
        		solve2(board,record,board.length-1,i);
        }
        for(int i=0;i<board.length;i++)
        	for(int j=0;j<board[0].length;j++){
        		if(board[i][j] == 'O' && record[i][j] == false)
        			board[i][j] = 'X';
        	}
    }

    //build binary tree from sorted linked list
    public TreeNode buildBinrayTree(ListNode head, ListNode tail){
    	if(head == null)	return null;
    	if(head == tail){
    		TreeNode n = new TreeNode(head.val);
    		return n;
    	}
    	int num = 1;
    	ListNode cur = head;
    	while(cur != tail){
    		cur = cur.next;
    		num ++;
    	}
    	int cnt = 0;
    	cur = head;
    	ListNode prev = null;
    	while(cnt < num/2){
    		prev = cur;
    		cur = cur.next;
    		cnt ++;
    	}
    	TreeNode left = buildBinrayTree(head,prev);
    	TreeNode right = null;
    	if(cur == tail)
    		right = null;
    	else
    		right = buildBinrayTree(cur.next,tail);
    	TreeNode n = new TreeNode(cur.val);
    	n.left = left;
    	n.right = right;
    	return n;
    }
    public TreeNode sortedListToBST(ListNode head) {
    	if(head == null)	return null;
        ListNode tail = head;
        while(tail.next != null)
        	tail = tail.next;
        return buildBinrayTree(head,tail);
    }

    //word ladder
    public int ladderLength(String start, String end, HashSet<String> dict) {
        if(start.equals(end))	return 1;
        if(dict.size() == 0)	return 0;
        Queue<String> q = new LinkedList<String>();
        q.offer(start);
        int prev_num = 1, cur_num = 0;
        int step = 1;
        while(!q.isEmpty()){
        	String cur = q.poll();
        	for(int i=0;i<cur.length();i++){
        		for(char a='a'; a<='z'; a++){
        			String new_str = cur.substring(0,i) + a + cur.substring(i+1,cur.length());
        			if(new_str.equals(end))		return step+1;
        			if(dict.contains(new_str)){
        				dict.remove(new_str);
        				q.offer(new_str);
        				cur_num ++;
        			}
        		}
        	}
        	prev_num --;
        	if(prev_num == 0){
        		step ++;
        		prev_num = cur_num;
        		cur_num = 0;
        	}
        }
        return 0;
    }

    //word ladder 2
    //BFS+STORE COMPLETE PATH  TLE
    //BFS+PARENT TABLE  TLE
    public void find_all_paths(HashMap<String,HashSet<String>> prev_dict, ArrayList<ArrayList<String>> result, String end, int step, String []record){
    	if(step == -1){
    		ArrayList<String> tmp = new ArrayList<String>();
    		for(int i=0;i<record.length;i++)
    			tmp.add(record[i]);
    		result.add(tmp);
    		return;
    	}
    	record[step] = end;
    	HashSet<String> prev_list = prev_dict.get(end);
    	for(String prevS : prev_list){
    		find_all_paths(prev_dict,result,prevS,step-1,record);
    	}
    }
    public ArrayList<ArrayList<String>> findLadders(String start, String end, HashSet<String> dict) {
        ArrayList<ArrayList<String>> result = new ArrayList<ArrayList<String>>();
        if(dict.size() == 0)	return result;
        if(start.equals(end))	{ ArrayList<String> tmp = new ArrayList<String>(); tmp.add(end); result.add(tmp); return result; }
        Queue<String> q = new LinkedList<String>();
        ArrayList<String> should_del = new ArrayList<String>();
        HashMap<String, HashSet<String>> prev_dict = new HashMap<String, HashSet<String>>();
        int prev_num = 1, cur_num = 0, step = 1;
        boolean isFound = false;
        q.offer(start);
        while(!q.isEmpty()){
        	String cur = q.poll();
        	prev_num --;
        	for(int i=0;i<cur.length();i++){
        		for(char c = 'a'; c <= 'z'; c++){
        			String new_str = cur.substring(0,i) + c + cur.substring(i+1,cur.length());
        			if(new_str.equals(end)){
        				 isFound = true;
        				 if(!prev_dict.containsKey(new_str)){
        					 HashSet<String> hs = new HashSet<String>();
        					 hs.add(cur);
        					 prev_dict.put(new_str, hs);
        				 }else	prev_dict.get(new_str).add(cur);
        			}
        			if(dict.contains(new_str)){
        				should_del.add(new_str);
        				q.offer(new_str);
        				cur_num ++;
        				if(!prev_dict.containsKey(new_str)){
        					HashSet<String> hs = new HashSet<String>();
        					hs.add(cur);
        					prev_dict.put(new_str, hs);
        				}else{
        					prev_dict.get(new_str).add(cur);
        				}
        			}
        		}
        	}
        	if(prev_num == 0){
        		prev_num = cur_num;
        		cur_num = 0;
        		step ++;
        		for(int i=0;i<should_del.size();i++){
        			if(dict.contains(should_del.get(i)))
        				dict.remove(should_del.get(i));
        		}
        		should_del.clear();
        		if(isFound){
        			String []record = new String[step];
        			find_all_paths(prev_dict,result,end,step-1,record);
        			return result;
        		}
        	}
        }
        return result;
    }
    
    //regular expression matching
    public boolean isMatch(String s, String p) {
        if(s.length() == 0 && p.length() == 0)		return true;
        if(s.length() == 0 && p.length() >= 2 && p.charAt(1) == '*')
        	return isMatch(s,p.substring(2,p.length()));
        if(s.length() == 0 || p.length() == 0)		return false;
        if(p.length() == 1){
        	char c = p.charAt(0);
        	if(c == '*')	return false;
        	if(c == '.' || c == s.charAt(0))	return isMatch(s.substring(1,s.length()),p.substring(1,p.length()));
        	return false;
        }else{
        	char c = p.charAt(0);
        	char next = p.charAt(1);
        	if(next == '*'){
        		if(c == '.'){
        			for(int i=0;i<=s.length();i++){
        				if(isMatch(s.substring(i,s.length()),p.substring(2,p.length())) == true)
        					return true;
        			}
        			return false;
        		}else{
        			if(isMatch(s,p.substring(2,p.length())) == true)	return true;
        			for(int i=0;i<s.length();i++){
        				if(s.charAt(i) == c){
        					if(isMatch(s.substring(i+1,s.length()),p.substring(2,p.length())) == true)
        						return true;
        				}else	break;
        			}
        			return false;
        		}
        	}else{
        		if(c == '.' || c == s.charAt(0))
        			return isMatch(s.substring(1,s.length()),p.substring(1,p.length()));
        		else
        			return false;
        	}
        }
    }

    //scramble string
    // think about how to write dp
    public boolean isScramble(String s1, String s2) {
    	if(s1.length() == 1)
    		return s1.charAt(0) == s2.charAt(0);
    	HashMap<Character,Integer> nums = new HashMap<Character,Integer>();
    	for(int i=0;i<s1.length();i++){
    		char c = s1.charAt(i);
    		if(!nums.containsKey(c))
    			nums.put(c, 1);
    		else
    			nums.put(c, nums.get(c)+1);
    	}
    	for(int i=0;i<s2.length();i++){
    		char c = s2.charAt(i);
    		if(!nums.containsKey(c) || nums.get(c) == 0)	return false;
    		nums.put(c, nums.get(c)-1);
    	}
    	for(int i=1;i<s1.length();i++){
    		//divide
    		String s1left = s1.substring(0,i);
    		String s1right = s1.substring(i,s1.length());
    		String s2left = s2.substring(0,i);
    		String s2right = s2.substring(i,s2.length());
    		boolean tmp_res = isScramble(s1left,s2left) && isScramble(s1right,s2right);
    		if(tmp_res == true)		return true;
    		else{
    			s2left = s2.substring(0,s2.length()-i);
    			s2right = s2.substring(s2.length()-i,s2.length());
    			tmp_res = isScramble(s1left,s2right) && isScramble(s1right,s2left);
    			if(tmp_res == true)		return true;
    		}
    	}
    	return false;
    }

    //wildcard matching
    //redo. Recur TLE
    public boolean isMatch2(String s, String p) {
    	boolean isStar = false;
    	int i=0,j=0;
    	int prevs = 0, prevp = 0;
    	while(i<s.length()){
    		if(j == p.length()){
    			if(!isStar)		return false;
    			else{
    				i = prevs + 1;
    				prevs ++;
    				j = prevp;
    			}
    		}else if(p.charAt(j) == '?'){
    			i ++;
    			j ++;
    		}else if(p.charAt(j) == '*'){
    			isStar = true;
    			prevs = i;
    			for(;j<p.length();j++){
    				if(p.charAt(j) != '*')	break;
    			}
    			if(j == p.length())		return true;
    			prevp = j;
    		}else{
    			if(s.charAt(i) != p.charAt(j)){
    				if(!isStar)
    					return false;
    				else{
    					i = prevs + 1;
    					prevs ++;
    					j = prevp;
    				}
    			}else{
    				i ++;
    				j ++;
    			}
    		}
    	}
    	for(;j<p.length();j++){
    		if(p.charAt(j) != '*')
    			return false;
    	}
    	return true;
    }
    
    //word search
    public boolean wordSearch(char [][]board, boolean [][]record, String word,int x,int y){
    	if(word.length() == 0)		return true;
    	if(x < 0 || x >= board.length || y < 0 || y >= board[0].length || record[x][y] || board[x][y] != word.charAt(0))		return false;
    	record[x][y] = true;
    	for(int i=-1;i<=1;i+=2){
    		if(wordSearch(board,record,word.substring(1,word.length()),x+i,y) == true)
    			return true;
    		if(wordSearch(board,record,word.substring(1,word.length()),x,y+i) == true)
    			return true;
    	}
    	record[x][y] = false;
    	return false;
    }
    public boolean exist(char[][] board, String word) {
        for(int i=0;i<board.length;i++){
        	for(int j=0;j<board[0].length;j++){
        		if(board[i][j] == word.charAt(0)){
        			boolean [][]record = new boolean[board.length][board[0].length];
        			if(wordSearch(board,record,word,i,j) == true)
        				return true;
        		}
        	}
        }
        return false;
    }

    //zigzag conversion
    public String convert(String s, int nRows) {
        String result = "";
        if(nRows == 1)		return s;
        char [][]tmp = new char[nRows][s.length()];
        int rowI = 0, colI = 0;
        boolean isDown = true;
        for(int i=0;i<s.length();i++){
        	char c = s.charAt(i);
        	if(isDown){
        		tmp[rowI++][colI] = c;
        		if(rowI == nRows){
        			isDown = false;
        			rowI -= 2;
        			colI ++;
        		}
        	}else{
        		tmp[rowI][colI] = c;
        		if(rowI == 0){
        			isDown = true;
        			rowI ++;
        		}else{
        			rowI --;
        			colI ++;
        		}
        	}
        }
        for(int i=0;i<nRows;i++){
        	for(int j=0;j<tmp[0].length;j++){
        		if(tmp[i][j] != 0)
        			result = result + tmp[i][j];
        	}
        }
        return result;
    }

    //integer to roman 
    public String intToRoman(int num) {
        String result = "";
        String one = "", five = "", ten = "";
        int bit_num = 0;
        while(num > 0){
        	int tmp = num % 10;
        	num = num / 10;
        	bit_num ++;
        	switch(bit_num){
        	case 1:
        		one = "I";
        		five = "V";
        		ten = "X";
        		break;
        	case 2:
        		one = "X";
        		five = "L";
        		ten = "C";
        		break;
        	case 3:
        		one = "C";
        		five = "D";
        		ten = "M";
        		break;
        	case 4:
        		one = "M";
        	}
        	if(tmp <= 3){
        		for(int i=1;i<=tmp;i++)
        			result = one + result;
        	}else if(tmp == 4){
        		result = one + five + result;
        	}else if(tmp <= 8){
        		for(int i=6;i<=tmp;i++)
        			result = one + result;
        		result = five + result;
        	}else{
        		result = one + ten + result;
        	}
        }
        return result;
    }

    //roman to integer
    public int romanToInt(String s) {
        int result = 0, i = 0;
        while(i != s.length()){
        	char c = s.charAt(i);
        	if(c == 'I'){
        		if(i+1 == s.length() || (s.charAt(i+1) != 'V' && s.charAt(i+1) != 'X')){
        			result = result + 1;
        			i ++;
        		}else if(s.charAt(i+1) == 'V'){
        			result = result + 4;
        			i += 2;
        		}else if(s.charAt(i+1) == 'X'){
        			result += 9;
        			i += 2;
        		}
        	}
        	else if(c == 'V'){
        		result += 5;
        		i ++;
        	}else if(c == 'X'){
        		if(i+1 == s.length() || (s.charAt(i+1) != 'L' && s.charAt(i+1) != 'C')){
        			result += 10;
        			i ++;
        		}else if(s.charAt(i+1) == 'L'){
        			result += 40;
        			i += 2;
        		}else if(s.charAt(i+1) == 'C'){
        			result += 90;
        			i += 2;
        		}
        	}else if(c == 'L'){
        		result += 50;
        		i ++;
        	}else if(c == 'C'){
        		if(i+1 == s.length() || (s.charAt(i+1) != 'D' && s.charAt(i+1) != 'M')){
        			result += 100;
        			i ++;
        		}else if(s.charAt(i+1) == 'D'){
        			result += 400;
        			i += 2;
        		}else if(s.charAt(i+1) == 'M'){
        			result += 900;
        			i += 2;
        		}
        	}else if(c == 'D'){
        		result += 500;
        		i ++;
        	}else if(c == 'M'){
        		result += 1000;
        		i++;
        	}
        }
        return result;
    }

    //valid sudoku
    public boolean isValidSudoku(char[][] board) {
        if(board.length %3 != 0 || board[0].length %3 != 0 || board.length == 0 || board[0].length == 0 || board.length != board[0].length)
        	return false;
        boolean []flag = new boolean[10];
        for(int i=0;i<board.length;i++){
        	for(int j=1;j<=9;j++)  flag[j] = false;
        	for(int j=0;j<board[0].length;j++){
        		if(board[i][j] == '.')		continue;
        		if(flag[board[i][j]-'0'] == true)
        			return false;
        		flag[board[i][j]-'0'] = true;
        	}
        }
        for(int i=0;i<board[0].length;i++){
        	for(int j=1;j<=9;j++)	flag[j] = false;
        	for(int j=0;j<board.length;j++){
        		if(board[j][i] == '.')		continue;
        		if(flag[board[j][i]-'0'] == true)	return false;
        		flag[board[j][i]-'0'] = true;
        	}
        }
        for(int i=0;i<board.length;i+=3){
        	for(int j=0;j<board[0].length;j+=3){
        		for(int k=1;k<=9;k++)		flag[k] = false;
        		for(int k1=0;k1<=2;k1++){
        			for(int k2=0;k2<=2;k2++){
        				if(board[i+k1][j+k2] == '.')	continue;
        				if(flag[board[i+k1][j+k2]-'0'] == true)		return false;
        				flag[board[i+k1][j+k2]-'0'] = true;
        			}
        		}
        	}
        }
        return true;
    }

    //sudoku solver
    public boolean solveSudokuRec(char [][]board, int i, int j, ArrayList<boolean []> row, ArrayList<boolean []>col, ArrayList<boolean []> grid){
    	if(i == board.length)	return true;
    	if(j == board[0].length)	return solveSudokuRec(board,i+1,0,row,col,grid);
    	if(board[i][j] != '.')	return solveSudokuRec(board,i,j+1,row,col,grid);
    	for(int num = 1; num <= 9; num++){
    		int grid_index = (i/3)*3+j/3;
    		if(row.get(i)[num] == true || col.get(j)[num] == true || grid.get(grid_index)[num] == true)
    			continue;
    		row.get(i)[num] = col.get(j)[num] = grid.get(grid_index)[num] = true;
    		board[i][j] = (char)(num+'0');
    		boolean tmp = solveSudokuRec(board,i,j+1,row,col,grid);
    		if(tmp == true)		return true;
    		board[i][j] = '.';
    		row.get(i)[num] = col.get(j)[num] = grid.get(grid_index)[num] = false;
    	}
    	return false;
    }
    public void solveSudoku(char[][] board) {
    	ArrayList<boolean []> row = new ArrayList<boolean []>();
    	ArrayList<boolean []> col = new ArrayList<boolean []>();
    	ArrayList<boolean []> grid = new ArrayList<boolean []>();
        for(int i=0;i<board.length;i++){
        	boolean []tmp = new boolean[10];
        	for(int j=0;j<board[0].length;j++){
        		if(board[i][j] == '.')		continue;
        		tmp[board[i][j]-'0'] = true;
        	}
        	row.add(tmp);
        }
        for(int i=0;i<board[0].length;i++){
        	boolean []tmp = new boolean[10];
        	for(int j=0;j<board.length;j++){
        		if(board[j][i] == '.')		continue;
        		tmp[board[j][i] - '0'] = true;
        	}
        	col.add(tmp);
        }
        for(int i=0;i<board.length;i+=3){
        	for(int j=0;j<board[0].length;j+=3){
        		boolean []tmp = new boolean[10];
        		for(int k1=0;k1<=2;k1++){
        			for(int k2=0;k2<=2;k2++){
        				if(board[i+k1][j+k2] == '.')	continue;
        				tmp[board[i+k1][j+k2] - '0'] = true;
        			}
        		}
        		grid.add(tmp);
        	}
        }
        solveSudokuRec(board,0,0,row,col,grid);
    }
    
    //multiply strings
    public String str_add(String num1, String num2){
    	String result = "";
    	int carry = 0;
    	if(num1.length() < num2.length()){
    		String tmp = num2;
    		num2 = num1;
    		num1 = tmp;
    	}
    	int j = num1.length()-1;
    	for(int i=num2.length()-1;i>=0;i--){
    		int tmp = num1.charAt(j) - '0' + num2.charAt(i) - '0' + carry;
    		j --;
    		int tmp2 = tmp % 10;
    		carry = tmp / 10;
    		result = tmp2 + result;
    	}
    	for(;j>=0;j--){
    		int tmp = num1.charAt(j) - '0' + carry;
    		int tmp2 = tmp % 10;
    		carry = tmp / 10;
    		result = tmp2 + result;
    	}
    	if(carry != 0)
    		result = carry + result;
    	return result;
    }
    public String mul_bit(String num1, int multiplied){
    	int carry = 0;
    	String result = "";
    	for(int i=num1.length()-1;i>=0;i--){
    		int a = num1.charAt(i) - '0';
    		int tmp = (a * multiplied + carry);
    		int tmp2 = tmp % 10;
    		carry = tmp / 10;
    		result = tmp2 + result;
    	}
    	if(carry != 0)
    		result = carry + result;
    	return result;
    }
    public String multiply(String num1, String num2) {
        String result = "0";
        int cur = 0;
        for(int i=num2.length()-1;i>=0;i--){
        	int multiplied = num2.charAt(i) - '0';
        	String tmp_result = mul_bit(num1,multiplied);
        	for(int j = 0; j < cur; j++)
        		tmp_result = tmp_result + "0";
        	result = str_add(result,tmp_result);
        	cur ++;
        }
        int i = 0;
        for(i=0;i<result.length();i++){
        	if(result.charAt(i) != '0')
        		return result.substring(i,result.length());
        }
        return "0";
    }
    
    //median of two sorted arrays
    public int findKthEle(int A[], int B[], int sa, int ea, int sb, int eb, int k){
    	if(sa > ea)		return B[sb+k-1];
    	if(sb > eb)		return A[sa+k-1];
    	if(k == 1)		return A[sa] < B[sb] ? A[sa] : B[sb];
    	int ma = (sa+ea)/2;
    	int mb = (sb+eb)/2;
    	int leftNum = ma - sa + mb - sb;	
    	if(A[ma] <= B[mb]){
    		if(leftNum + 1 >= k){
    			return findKthEle(A,B,sa,ea,sb,mb-1,k);
    		}else{
    			return findKthEle(A,B,ma+1,ea,sb,eb,k-(ma-sa)-1);
    		}
    	}else{
    		if(leftNum + 1 >= k){
    			return findKthEle(A,B,sa,ma-1,sb,eb,k);
    		}else{
    			return findKthEle(A,B,sa,ea,mb+1,eb,k-(mb-sb)-1);
    		}
    	}
    }
    public double findMedianSortedArrays(int A[], int B[]) {
        if((A.length+B.length)%2 == 1)
        	return findKthEle(A,B,0,A.length-1,0,B.length-1,(A.length+B.length)/2+1);
        else
        	return ((double)findKthEle(A,B,0,A.length-1,0,B.length-1,(A.length+B.length)/2) + (double)findKthEle(A,B,0,A.length-1,0,B.length-1,(A.length+B.length)/2+1))/2;
    }
    
    public static void main(String args[]){
    	int []A = {1,3,5,8};
    	int []B = {2,4,6,7};
    	System.out.println(new practice().findMedianSortedArrays(A, B));
	}
}
