package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class CoreTest {

  private Core mCore;
  
  @Before
  public void setUp() throws Exception {
    // init logging
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    // create core
    mCore = new Core();
  }
  
  @Test
  public void testLastUpdateCount(){
    
    // create graph with 10 vertices
    final DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    for (int i=0; i<10; i++) graph.addVertex(new ScalarVertex(i));
    
    // schedule a simple updater on all 10 vertices
    mCore.setGraph(graph);
    mCore.scheduleAll(new DefaultUpdater(){
      @Override
      public void update(Context context, ScalarVertex vertex){}

      @Override
      protected DefaultUpdater clone() {
        return this;
      }
    });
    mCore.start();
    
    // check count
    assertEquals("Checking lastUpdateCount", 10, mCore.lastUpdateCount());
    
  }
  
  @Test
  public void testGlobalConst(){
    
    // set const
    Integer obj = 6;
    mCore.addGlobalConst("obj", obj);
    
    // do some computation
    // create graph with 10 vertices
    final DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    for (int i=0; i<10; i++) graph.addVertex(new ScalarVertex(i));
    
    // schedule a simple updater on all 10 vertices
    mCore.setGraph(graph);
    mCore.scheduleAll(new DefaultUpdater(){
      @Override
      public void update(Context context, ScalarVertex vertex){
        // get const
        Integer returned = mCore.getGlobal("obj", Integer.class);
        assertEquals("Checking value of stored constant.", new Integer(6), returned);
      }

      @Override
      protected DefaultUpdater clone() {
        return this;
      }
    });
    mCore.start();
    
    // get const
    Integer returned = mCore.getGlobal("obj", Integer.class);
    assertEquals("Checking value of stored constant.", new Integer(6), returned);
    
  }
  
  @Test
  public void testGlobal(){
    
    // set const
    Integer obj = 6;
    mCore.addGlobal("obj", obj);
    
    // do some computation
    // create graph with 10 vertices
    final DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    for (int i=0; i<10; i++) graph.addVertex(new ScalarVertex(i));
    
    // schedule a simple updater on all 10 vertices
    mCore.setGraph(graph);
    mCore.scheduleAll(new DefaultUpdater(){
      @Override
      public void update(Context context, ScalarVertex vertex){
        // set value to 7
        mCore.setGlobal("obj", new Integer(7));
      }

      @Override
      protected DefaultUpdater clone() {
        return this;
      }
    });
    mCore.start();
    
    // get const
    Integer returned = mCore.getGlobal("obj", Integer.class);
    assertEquals("Checking value of stored constant.", new Integer(7), returned);
    
  }

  @After
  public void tearDown() throws Exception {
    mCore.destroy();
  }
  
  private abstract class DefaultUpdater extends Updater<ScalarVertex, DefaultWeightedEdge, DefaultUpdater> {}

}
