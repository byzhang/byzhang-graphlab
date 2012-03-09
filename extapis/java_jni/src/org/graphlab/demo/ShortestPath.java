package org.graphlab.demo;

import java.io.IOException;
import java.util.Scanner;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Aggregator;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Core.CoreException;
import org.graphlab.CoreConfiguration;
import org.graphlab.Scheduler;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.graphlab.util.GraphLoader;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Shortest Path algorithm.
 * 
 * <p>
 * Demonstrates GraphLab over JNI. To run this class, provide a path to a tsv
 * file in the program arguments. Some examples are available in the
 * <tt>java_jni/test-graphs</tt> directory.
 * </p>
 * 
 * <p>
 * For example:
 * 
 * <pre>
 * cd extapis/java_jni
 * java \
 *  -classpath bin:lib/log4j-1.2.16.jar \
 *  -Djava.library.path=../../release/src/graphlab/jni \
 *  org.graphlab.demo.ShortestPath test-graphs/toy.tsv
 * </pre>
 * 
 * </p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class ShortestPath {

  private static final Logger logger = Logger.getLogger(ShortestPath.class);

  public static void main(String[] args) {

    initLogger();
    
    Scanner sc = new Scanner(System.in);
    sc.next();
    
    logger.trace("Main method in " +
        ShortestPath.class.getCanonicalName() +
        " started.");

    // check arguments
    if (!checkParams(args)) {
      logger.trace("Exiting main method.");
      return;
    }

    String filename = args[0];
    logger.info("Graph file: " + filename);

    // initialize graphlab core
    final Core core;
    try {
      logger.trace("Initializing GraphLab core ...");
      CoreConfiguration config = new CoreConfiguration();
      config.setScheduler(Scheduler.FIFO);
      core = new Core();
    } catch (CoreException e) {
      logger.fatal("Unable to initialize core. Terminating.", e);
      logger.trace("Exiting main method.");
      return;
    }

    // construct graph
    logger.trace("Constructing graph from " + filename + " ...");
    final DirectedGraph<ScalarVertex, DefaultWeightedEdge> graph;
    try {
      graph = constructGraph(filename);
    } catch (IOException e) {
      logger.fatal("Unable to construct graph. Terminating.", e);
      core.destroy();
      logger.trace("Exiting main method.");
      return;
    }

    for (ScalarVertex v : graph.vertexSet()) {
      v.setValue(Integer.MAX_VALUE);
    }

    // start from root
    ScalarVertex root = graph.vertexSet().iterator().next();
    root.setValue(0);

    core.setGraph(graph);
    core.schedule(root, new ShortestPathUpdater(graph));

    logger.trace("Running graphlab ...");
    logger.trace("Runtime: " + core.start() + " s");
    logger.trace("Update count: " + core.lastUpdateCount());
    

    core.addAggregator("agg", new ShortestPathAggregator(), 0);
    core.aggregateNow("agg");
    
    // destroy core
    logger.trace("Destroying core ...");
    core.destroy();

  }

  private static void initLogger() {
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.ALL);
    logger.setLevel(Level.ALL);
  }

  /**
   * Checks that required input parameters are available and valid. Prints
   * instructions if not all parameters were valid.
   * 
   * @param args
   *          array of program arguments
   * @return true if parameters are OK; false otherwise
   */
  private static boolean checkParams(String[] args) {

    if (args.length != 1) {
      System.out.println("Please provide filename.");
      System.out.println("Usage: java -Djava.library.path=... "
          + ShortestPath.class.getCanonicalName() + " path/to/tsv/file");
      return false;
    }

    return true;

  }

  private static DirectedGraph<ScalarVertex, DefaultWeightedEdge>
    constructGraph(String filename)
    throws IOException {

    DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph =
       new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    GraphLoader.loadGraphFromTsvFile(graph, filename);
    return graph;

  }
  
  /**
   * Finds shortest path to a particular vertex
   * @author Jiunn Haur Lim
   */
  private static class ShortestPathUpdater extends Updater<ScalarVertex, DefaultWeightedEdge, ShortestPathUpdater> {

    private DirectedGraph<ScalarVertex, DefaultWeightedEdge> mGraph;

    public ShortestPathUpdater(DirectedGraph<ScalarVertex, DefaultWeightedEdge> g) {
      this.mGraph = g;
    }

    @Override
    public void update(Context context, ScalarVertex vertex) {
      
      // find shortest known distance into this node
      for (DefaultWeightedEdge edge : mGraph.incomingEdgesOf(vertex)) {
        ScalarVertex neighbor = mGraph.getEdgeSource(edge);
        double weight = mGraph.getEdgeWeight(edge);
        vertex.setValue(
            (float) Math.min(vertex.value(), neighbor.value() + weight)
        );
      }

      // reschedule any affected neighbors
      for (DefaultWeightedEdge edge : mGraph.outgoingEdgesOf(vertex)) {
        ScalarVertex neighbor = mGraph.getEdgeTarget(edge);
        double weight = mGraph.getEdgeWeight(edge);
        if (neighbor.value() > (vertex.value() + weight)) {
          context.schedule(neighbor, this);
        }
      }

    }

    @Override
    protected ShortestPathUpdater clone() {
      return new ShortestPathUpdater(mGraph);
    }

  }

  /**
   * Aggregates shortest path information: total distance, longest distance, and furthest vertex
   * @author Jiunn Haur Lim
   */
  private static class ShortestPathAggregator extends Aggregator<ScalarVertex, ShortestPathAggregator> {
  
    /** maximum distance */
    private double mMaxDist;
    
    /** total distance */
    private double mSumDist;
    
    /** furthest vertex */
    private int mFurthestVertex;
    
    @Override
    protected void exec(Context context, ScalarVertex vertex) {
      double dist = vertex.value();
      mSumDist += dist;
      if (dist > mMaxDist){
        mMaxDist = dist;
        mFurthestVertex = vertex.id();
      }
    }
    
    @Override
    protected void add (ShortestPathAggregator other) {
      // add distances
      mSumDist += other.mSumDist;
      // take max
      if (other.mMaxDist > mMaxDist){
        mMaxDist = other.mMaxDist;
        mFurthestVertex = other.mFurthestVertex;
      }
    }
    
    @Override
    protected void finalize(Context context) {
       // output results
       logger.info ("Total Distance:\t\t" + mSumDist + "\n" + 
                    "Longest Distance:\t" + mMaxDist + "\n" +
                    "Furthest Vertex:\t"  + mFurthestVertex + "\n");
    }

    @Override
    protected ShortestPathAggregator clone() {
      return new ShortestPathAggregator();
    }
    
  }
  
}
