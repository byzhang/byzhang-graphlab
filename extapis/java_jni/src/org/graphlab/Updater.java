package org.graphlab;

import org.graphlab.data.Vertex;

/**
 * Updater
 * 
 * <p>
 * The GraphLab engine will invoke an updater on each scheduled node. Extend
 * this class to provide an update function for that node. Note that the update
 * function may update node data, modify edge data, and schedule neighbors, but
 * may not modify the graph structure. You may reuse the updater object on
 * across multiple vertices (this is encouraged).
 * </p>
 * 
 * @param <V>
 *          Vertex type that will be used in {@link #update(Context, Vertex)}
 * @param <E>
 *          Edge type
 * @param <U>
 *          Updater
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public abstract class Updater<V extends Vertex, E, // TODO: easier generics?
U extends Updater<V, E, U>> {

  static {
    initNative();
  }

  /*
   * The set of edges that are operated on during gather and scatter
   */
  protected static final int IN_EDGES = 0;
  protected static final int OUT_EDGES = 1;
  protected static final int ALL_EDGES = 2;
  protected static final int NO_EDGES = 3;

  /** No locks are acquired at all. */
  protected static final int NULL_CONSISTENCY = 0;

  /** Write to self. No lock on adjacent. */
  protected static final int VERTEX_CONSISTENCY = 1;

  /** Write to self. Read from adjacent structures. */
  protected static final int EDGE_CONSISTENCY = 2;

  /** Write to self and adjacent structures. */
  protected static final int FULL_CONSISTENCY = 3;

  /** Use externally described default. */
  protected static final int DEFAULT_CONSISTENCY = 4;

  /**
   * Updates the vertex. Subclasses may wish to maintain a reference to the
   * graph object.
   * 
   * @param context
   *          graphlab context; use {@link Context#schedule(Vertex, Updater)} to
   *          schedule vertices.
   * @param vertex
   *          vertex to be updated
   */
  protected void update(Context context, V vertex) {
  }

  /**
   * When multiple update functors are scheduled to be run on the same function
   * they are added. The default behavior is to simply ignore the later update
   * functors. Override this method to implement your own behavior.
   * 
   * @param updater
   */
  protected void add(U updater) {
    return;
  }

  /**
   * Get the priority of the update functor. Defaults to 0.
   * 
   * @return priority
   */
  protected double priority() {
    return 0;
  }

  /**
   * Returns true if the factorized (gather, apply, scatter) version of the
   * update functor is to be used. False by default, so override this!
   * 
   * @return true if the factorized version is to be used; false otherwise
   */
  protected boolean isFactorizable() {
    return false;
  }

  /**
   * Override to return one of {@link #IN_EDGES}, {@link #OUT_EDGES},
   * {@link #ALL_EDGES}, and {@link #NO_EDGES}.
   * 
   * @return the set of edges to gather
   */
  protected int gatherEdges() {
    return IN_EDGES;
  }

  /**
   * Override to return one of {@link #IN_EDGES}, {@link #OUT_EDGES},
   * {@link #ALL_EDGES}, and {@link #NO_EDGES}.
   * 
   * @return the set of edges to scatter
   */
  protected int scatterEdges() {
    return OUT_EDGES;
  }

  /**
   * Gets the context range required by this update functor. If not implemented
   * by the derived class then the default context range is returned.
   * 
   * @return one of {@link #VERTEX_CONSISTENCY}, {@link #EDGE_CONSISTENCY},
   *         {@link #FULL_CONSISTENCY}, {@link #NULL_CONSISTENCY}, and
   *         {@link #DEFAULT_CONSISTENCY}
   */
  protected int consistency() {
    return DEFAULT_CONSISTENCY;
  }

  /**
   * Returns consistency required during gather ops.
   * 
   * @return one of {@link #VERTEX_CONSISTENCY}, {@link #EDGE_CONSISTENCY},
   *         {@link #FULL_CONSISTENCY}, {@link #NULL_CONSISTENCY}, and
   *         {@link #DEFAULT_CONSISTENCY}
   */
  protected int gatherConsistency() {
    return DEFAULT_CONSISTENCY;
  }

  /**
   * Returns consistency required during scatter ops.
   * 
   * @return one of {@link #VERTEX_CONSISTENCY}, {@link #EDGE_CONSISTENCY},
   *         {@link #FULL_CONSISTENCY}, {@link #NULL_CONSISTENCY}, and
   *         {@link #DEFAULT_CONSISTENCY}
   */
  protected int scatterConsistency() {
    return DEFAULT_CONSISTENCY;
  }

  /**
   * Init gather is called before gathering.
   */
  protected void initGather() {
    throw new UnsupportedOperationException();
  }

  /**
   * Gather is called on all gather_edges() and may be called in parallel. The
   * merge() operation is used to join update functors.
   * 
   * @param edge
   */
  protected void gather(E edge) {
    throw new UnsupportedOperationException();
  }

  /**
   * Merges update functors during the gather process.
   * 
   * @param updater
   *          updater from another gather operation
   */
  protected void merge(U updater) {
    throw new UnsupportedOperationException();
  }

  /**
   * Apply is called within the vertex consistency model on the center vertex
   * after all gathers have completed.
   * 
   * @param vertex
   *          vertex to operate on
   */
  protected void apply(V vertex) {
    throw new UnsupportedOperationException();
  }

  /**
   * Scatter is invoked on all scatter_edges() after calling and may be invoked
   * in parallel.
   * 
   * @param context
   *          use this to schedule neighbors
   * @param edge
   */
  protected void scatter(Context context, E edge) {
    throw new UnsupportedOperationException();
  }

  /*
   * Required because multiple updaters might be executed in parallel.
   * (non-Javadoc)
   * 
   * @see java.lang.Object#clone()
   */
  @Override
  protected abstract U clone();

  /**
   * Executes the updater on the specified vertex. This is <em>only</em> invoked
   * by the proxy updater in the JNI library.
   * 
   * @param contextPtr
   *          address of graphlab::icontext_type object
   * @param vertexId
   *          application vertex ID
   */
  private void update(long contextPtr, V vertex) {
    Context context = new Context(contextPtr);
    update(context, vertex);

  }

  /**
   * Executes scatter on the specified edge. This is <em>only</em> invoked by
   * the proxy updater in the JNI library.
   * 
   * @param contextPtr
   * @param edge
   */
  private void scatter(long contextPtr, E edge) {
    Context context = new Context(contextPtr);
    scatter(context, edge);
  }

  /**
   * Initialize native class (set field IDs and method IDs)
   */
  private static native void initNative();

}
