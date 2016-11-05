package edu.koritsas.slimemold.shortestpath;

import edu.koritsas.slimemold.AbstractDirectedPhysarumPolycephalum;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Node;

/**
 * Created by ilias on 12/10/2016.
 */
public abstract class DirecredPhysarumPolycephalumSP extends AbstractDirectedPhysarumPolycephalum {
    public DirecredPhysarumPolycephalumSP(DirectedGraph graph, Node sourceNode, Node sinkNode, int numberOfIterations) {
        super(graph, sourceNode, sinkNode, numberOfIterations);
    }
}
