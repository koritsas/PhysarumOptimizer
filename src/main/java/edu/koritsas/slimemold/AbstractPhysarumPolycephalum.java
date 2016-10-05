package edu.koritsas.slimemold;

import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.jfree.data.xy.XYSeries;

import java.util.HashMap;
import java.util.List;
import java.util.logging.Logger;

/**
 * Created by ilias on 5/10/2016.
 */
public abstract class AbstractPhysarumPolycephalum {
    private Graph graph;
    private Node sourceNode;
    private Node sinkNode;
    private List<Node> sinkNodesList;
    protected final double Io;
    public final double γ;
    protected final int numberOfIterations;
    protected Node sinkNode;
    protected Node sourceNode;
    private HashMap<Edge,Double> fluxMap ;
    private HashMap<Edge,Double> conductivityMap;
    private HashMap<Node,Double> pressureMap;

    protected static final Logger logger = Logger.getLogger(PhysarumPolycephalum.class.getName());

    protected int iteration=0;

    private HashMap<Edge, XYSeries> flowSeriesMap = new HashMap<Edge, XYSeries>();
    private HashMap<Edge, XYSeries> conductivitySeriesMap =new HashMap<Edge, XYSeries>();
    private HashMap<Node, XYSeries> pressureSeriesMap =new HashMap<Node, XYSeries>();

    public HashMap<Edge, XYSeries> getFlowSeriesMap() {
        return flowSeriesMap;
    }

    public void setFlowSeriesMap(HashMap<Edge, XYSeries> flowSeriesMap) {
        this.flowSeriesMap = flowSeriesMap;
    }

    public HashMap<Edge, XYSeries> getConductivitySeriesMap() {
        return conductivitySeriesMap;
    }

    public void setConductivitySeriesMap(HashMap<Edge, XYSeries> conductivitySeriesMap) {
        this.conductivitySeriesMap = conductivitySeriesMap;
    }




}
