package edu.koritsas.slimemold.shortestpath;

import edu.koritsas.slimemold.AbstractDirectedPhysarumPolycephalum;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

/**
 * Created by ilias on 12/10/2016.
 */
public abstract class DirecredPhysarumPolycephalumSP extends AbstractDirectedPhysarumPolycephalum {

    public DirecredPhysarumPolycephalumSP(Graph graph, Node sourceNode, Node sinkNode, double absoluteThreshold, double relativeThreshold, int numberOfIterations) {
        super(graph, sourceNode, sinkNode, absoluteThreshold, relativeThreshold, numberOfIterations);
    }
}
