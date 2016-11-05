package edu.koritsas.slimemold.shortestpath;

import edu.koritsas.slimemold.AbstractPhysarumPolycephalum;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

/**
 * Created by ilias on 5/11/2016.
 */
public  abstract  class PhysarumPolycephalumCSP extends AbstractPhysarumPolycephalum {
    public PhysarumPolycephalumCSP(Graph graph, Node sourceNode, Node sinkNode, int numberOfIterations) {
        super(graph, sourceNode, sinkNode, numberOfIterations);
    }
}
