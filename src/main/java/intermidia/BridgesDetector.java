/******************************************************************************
 *
 *  Identifies bridge edges and prints them out. This decomposes
 *  a directed graph into two-edge connected components.
 *  Runs in O(E + V) time.
 *
 *  Key quantity:  low[v] = minimum DFS preorder number of v
 *  and the set of vertices w for which there is a back edge (x, w)
 *  with x a descendant of v and w an ancestor of v.
 *
 *  Note: code assumes no parallel edges, e.g., two parallel edges
 *  would be (incorrectly) identified as bridges.
 *  
 *  code from: http://algs4.cs.princeton.edu/ 
 *  
 *  Note: adapted to work with GraphStream API
 *
 ******************************************************************************/

package intermidia;

import java.util.ArrayList;

import org.graphstream.graph.Edge;
import org.graphstream.graph.implementations.SingleGraph;
import org.openimaj.util.pair.IntIntPair;

public class BridgesDetector {
    private int bridgeQuantity;      // number of bridges
    private int cnt;          // counter
    private int[] pre;        // pre[v] = order in which dfs examines v
    private int[] low;        // low[v] = lowest preorder of any vertex connected to v
    private ArrayList<IntIntPair> bridges;

    public BridgesDetector(SingleGraph G) {
        low = new int[G.getNodeCount()];
        pre = new int[G.getNodeCount()];
        bridges = new ArrayList<IntIntPair>();
        for (int v = 0; v < G.getNodeCount(); v++)
            low[v] = -1;
        for (int v = 0; v < G.getNodeCount(); v++)
            pre[v] = -1;
        
        for (int v = 0; v < G.getNodeCount(); v++)
            if (pre[v] == -1)
                dfs(G, Integer.toString(v), Integer.toString(v));
    }

    public int componentQuantity() { return bridgeQuantity + 1; }
    
    public ArrayList<IntIntPair> bridges()
    {
    	return bridges;
    }

    private void dfs(SingleGraph G, String su, String sv) 
    {
    	int u = Integer.parseInt(su);
    	int v = Integer.parseInt(sv);
        pre[v] = cnt++;
        low[v] = pre[v];
        for (Edge e : G.getNode(sv).getEachEdge()) 
        {   
        	//Get ending node of the edge
        	String sw = e.getNode0().getId();
        	if(sw.compareTo(sv) == 0)
        	{
        		sw = e.getNode1().getId();
        	}
        	int w = Integer.parseInt(sw);
        	
            if (pre[w] == -1) 
            {
                dfs(G, sv, sw);
                low[v] = Math.min(low[v], low[w]);
                if (low[w] == pre[w]) 
                {
                	bridges.add(new IntIntPair(v, w));
                    bridgeQuantity++;
                }
            }

            // update low number - ignore reverse of edge leading to v
            else if (w != u)
                low[v] = Math.min(low[v], pre[w]);
        }
    }
}
