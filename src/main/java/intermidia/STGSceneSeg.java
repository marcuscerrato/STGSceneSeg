package intermidia;

import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.SingleGraph;
import org.openimaj.feature.DoubleFVComparison;
import org.openimaj.math.statistics.distribution.Histogram;
import org.openimaj.util.pair.IntIntPair;

import TVSSUnits.Shot;
import TVSSUnits.ShotList;

import com.opencsv.CSVReader;

public class STGSceneSeg 
{
	private static double delta = 0.15;
	private static int timeThreshold = 10;
		
	//Usage: <in: feature vectors> <out: segmentation result> <in: delta> <in: time window> <in: output format>
    public static void main( String[] args ) throws Exception
    {
    	//Read feature word histograms features from CSV file.
    	CSVReader featureReader = new CSVReader(new FileReader(args[0]), ' ');
		String [] line;
		ShotList shotList = new ShotList();
		int shotQty = 0;
		
		delta = Double.parseDouble(args[2]);
		timeThreshold = Integer.parseInt(args[3]);
		//System.out.println("Delta value: "+ delta + " / Time window value: " + timeThreshold);
		
		//Build shot list with feature words histograms
		while ((line = featureReader.readNext()) != null) 
		{
			shotList.addShot(new Shot());
			
			int fvSize = line.length - 1;
			double fv[] = new double[fvSize];
			
			double sumZero = 0;
			for(int i = 0; i < fvSize; i++)
			{
				//fv[i] = Integer.parseInt(line[i + 1]);
				fv[i] = Double.parseDouble(line[i + 1]);
				sumZero += fv[i];
			}
			
			//If the histogram sums zero put extreme values in it to avoid NaN when computing similarity
			if(sumZero == 0)
			{
				for(int i = 0; i < fvSize; i++)
				{
					fv[i] = Double.MIN_VALUE;
				}				
			}
			
			shotList.getLastShot().setFeatureWordHistogram(new Histogram( (new Histogram(fv)).normaliseFV() ));
			shotQty++;
		}
		featureReader.close();	
		
		
		double interShotDist[][] = new double[shotQty][shotQty];
		for(int i = 0; i < shotQty; i++)
		{
			for(int j = 0; j < shotQty; j++)
			{
				// 1 - for cosine distance
				interShotDist[i][j] =  1 - shotList.getShot(i).getFeatureWordHistogram().compare(shotList.getShot(j).getFeatureWordHistogram(), DoubleFVComparison.COSINE_SIM);
				
				/*For some unknown reason, Openimaj comparison sometimes returns 
				  negative and greater than 1 values.. so this code corrects it.*/
				if(interShotDist[i][j] < 0)
				{
					interShotDist[i][j] = 0;
				}
				if(interShotDist[i][j] > 1)
				{
					interShotDist[i][j] = 1;
				}
			}
		}
		
		ArrayList<ArrayList<Integer>> shotClusters = completeLinkThresholdClusterer(interShotDist, delta, timeThreshold);

	
		//STG Construction
		SingleGraph g = new SingleGraph("STG");
		for(int i = 0; i < shotClusters.size(); i++)
		{
			g.addNode(Integer.toString(i));
			g.getNode(Integer.toString(i)).setAttribute("label", i); //Visualization purposes only
		}
		
		//Combine all shot clusters
		for(int i = 0; i < shotClusters.size(); i++)
		{
			for(int j = i + 1; j < shotClusters.size(); j++)
			{
				/*Insert an edge from U to W (i and j for iteration purpose) if
				 	there is a Sl in U and a Sm in W such that m = l + 1. */  
				if(containsSubsequentShots(shotClusters.get(i), shotClusters.get(j)))
				{
					g.addEdge(Integer.toString(i) + "/" + Integer.toString(j), 
							Integer.toString(i), Integer.toString(j));
				}
			}
		}
		
		ConnectedComponents connectedComponents = new ConnectedComponents(g);
		connectedComponents.setCountAttribute("componentNum");

		BridgesDetector bridgesDetector = new BridgesDetector(g);
		
		ArrayList<IntIntPair> bridges = bridgesDetector.bridges();		
		for(IntIntPair bridge: bridges)
		{
			g.removeEdge(Integer.toString(bridge.first), Integer.toString(bridge.second));

		}	
		int componentQuantity = connectedComponents.getConnectedComponentsCount();
		
		//Divide shots into scenes
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] scenes = new ArrayList[componentQuantity];
		for(int i = 0; i < componentQuantity; i++)
		{
			scenes[i] = new ArrayList<Integer>();
		}
		
		Collection<Node> nodeSet = g.getNodeSet();
		for(Node n : nodeSet)
		{
			scenes[(int) n.getAttribute("componentNum")].add(Integer.parseInt(n.getId()));			
		}				
		
		ArrayList<SceneShots> sceneSeg = new ArrayList<SceneShots>();
		for(int i = 0; i < componentQuantity; i++)
		{			
			ArrayList<Integer> allSceneShots = new ArrayList<Integer>();
			for(Integer n : scenes[i])
			{
				allSceneShots.addAll(shotClusters.get(n));
			}
			SceneShots boundaries = findBoundaries(allSceneShots);			
			sceneSeg.add(boundaries);		
		}
		Collections.sort(sceneSeg);
		
		int lastAnalyzedShotIndex = 0;
		int lastWrittenShotIndex = 0;
    	FileWriter sceneWriter = new FileWriter(args[1]);
    	
    	for(SceneShots sceneShotsPair : sceneSeg)
    	{
    		/*Heuristic 1 to prevent scenes contained inside other scenes*/
    		if( !(sceneShotsPair.first < lastAnalyzedShotIndex) )
    		{
    			/*Heuristic 2 to prevent gaps between scenes*/
    			if(sceneShotsPair.first > (lastWrittenShotIndex + 1))
    			{
    				
        			/*Trojahn output standard*/
        			if(args.length > 4 && args[4].compareTo("-alt_csv_out") == 0 )
        			{
        				sceneWriter.write((lastWrittenShotIndex + 2) + ", " + (sceneShotsPair.first) + "\n");
        			}
        			/*Simple csv output standard*/
        			else
        			{
        				sceneWriter.write((lastWrittenShotIndex + 1) + " " + (sceneShotsPair.first - 1) + "\n");
        			}	
    			}
    			/*Heuristic 2 ends here*/
    			
    			
    			/*Trojahn output standard*/
    			if(args.length > 4 && args[4].compareTo("-alt_csv_out") == 0 )
    			{
    				sceneWriter.write((sceneShotsPair.first + 1) + ", " + (sceneShotsPair.last + 1) + "\n");
    			}
    			/*Simple csv output standard*/
    			else
    			{
    				sceneWriter.write(sceneShotsPair.first + " " + sceneShotsPair.last + "\n");
    			}
    			lastWrittenShotIndex = sceneShotsPair.last;
    		}
    		/*Part of Heuristic 1*/
    		if(lastAnalyzedShotIndex < sceneShotsPair.last)
    		{
    			lastAnalyzedShotIndex = sceneShotsPair.last;
    		}    		
    	}
    	/*Heuristic 3 prevents loosing final shots*/
    	if(lastWrittenShotIndex < shotQty - 1)
    	{
			/*Trojahn output standard*/
			if(args.length > 4 && args[4].compareTo("-alt_csv_out") == 0 )
			{
				sceneWriter.write((lastWrittenShotIndex + 2) + ", " + shotQty + "\n");
			}
			/*Simple csv output standard*/
			else
			{
				sceneWriter.write((lastWrittenShotIndex + 1) + " " + (shotQty - 1) + "\n");
			}	
    	}
    	/*Heuristic 3 ends here*/
    	
    	sceneWriter.close();
    	//System.out.println("Segmentation process ended.");
    	//System.exit(0);
		
		
    }
    
    //Allows comparation of IntIntPairs with the implementation of compareTo method
    public static class SceneShots implements Comparable<SceneShots>
    {
    	public int first;
    	public int last;
    	
    	public SceneShots(int first, int second) 
    	{
			this.first = first;
			this.last = second;
		}
    	
    	
		@Override
		public int compareTo(SceneShots o) 
		{
			return this.first - o.first;
		}
    	
    }
    
    
    //Verifies if there is a subsequent shot pair of shots in U, W.
    private static boolean containsSubsequentShots(ArrayList<Integer> u, ArrayList<Integer> w)
    {
    	boolean returnValue = false;
    	for(Integer sl : u)
    	{
    		for(Integer sm : w)
    		{
    			if( (sl == (sm + 1)) || (sm == (sl + 1)) )
    			{
    				return true;
    			}
    		}
    	}    	
    	return returnValue;
    }
    
    private static ArrayList<ArrayList<Integer>> completeLinkThresholdClusterer(double [][] distances, double delta, int timeThreshold)
    {
    	ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>();    	

    	//Initialize unitary clusters
    	for(int i = 0; i < distances.length; i++)
    	{
    		clusters.add(new ArrayList<Integer>());
    		clusters.get(i).add(i);
    	}
    	    	
    	//While there are clusters with distances greater than delta threshold and there are more than one cluster
    	while( !allClustersDistsGreaterThanDelta(clusters, delta, distances, timeThreshold) && (clusters.size() > 1))
    	{    	
    		//Find minimum distance between all clusters
    		int clusterQty = clusters.size();    		
    		double minDist = distances[0][1]; //Initialize with any valid distance
    		int clusterR = 0, clusterS = 1; //Must be compatible with the minDist initialization
    		for(int i = 0; i < clusterQty; i++)
    		{
    			for(int j = i + 1; j < clusterQty; j++)
    			{
    				double actualDist = clusterDist(clusters.get(i), clusters.get(j), distances, timeThreshold);
    				if(actualDist < minDist)
    				{
    					minDist = actualDist;
    					clusterR = i;
    					clusterS = j;
    				}
    			}
    		}
    		//Merge R and S into a new cluster
    		for(Integer shot: clusters.get(clusterS))
    		{
    			clusters.get(clusterR).add(shot);
    		}
    		clusters.remove(clusterS);
    	}
    		
    	return clusters;    	
    }
    
    //Returns the smaller and greater indexes (of shots) in the input list.
    private static SceneShots findBoundaries(ArrayList<Integer> shotList)
    {
    	if(shotList == null || shotList.isEmpty())
    		return null;
    	
    	int min, max;
    	Iterator<Integer> shotListIterator = shotList.iterator();
    	min = shotListIterator.next();
    	max = min;
    	while(shotListIterator.hasNext())
    	{
    		int next = shotListIterator.next();
    		if(next < min)
    		{
    			min = next;
    		}
    		if(next > max)
    		{
    			max = next;
    		}
    	}
    	return new SceneShots(min, max);
    }
    
    //Find maximum distance between all clusters
    private static boolean allClustersDistsGreaterThanDelta(ArrayList<ArrayList<Integer>> clusters, double delta, double[][] distances, int timeThreshold)
    {    	
		int clusterQty = clusters.size();
		for(int i = 0; i < clusterQty; i++)
		{
			for(int j = i + 1; j < clusterQty; j++)
			{
				if( clusterDist(clusters.get(i), clusters.get(j), distances, timeThreshold) <= delta )
				{				
					return false;
				}
			}
		}    		
    	return true;
    }
    
    //Compute distance between two clusters of shots
    private static double clusterDist(ArrayList<Integer> clusterA, ArrayList<Integer> clusterB, double[][] distances, int timeThreshold)
    {
    	double maxDist = 0;
    	for(Integer a : clusterA)
    	{
    		for(Integer b : clusterB)
    		{
    			if(distances[a][b] > maxDist)
    			{
    				//Applies the temporal distance criterion    				
    				if(Math.abs(a - b) <= timeThreshold)
    				{
    					maxDist = distances[a][b];
    				}else
    				{
    					maxDist = Double.MAX_VALUE;
    				}
    			}
    		}
    	}
    	return maxDist;
    }
                
}
