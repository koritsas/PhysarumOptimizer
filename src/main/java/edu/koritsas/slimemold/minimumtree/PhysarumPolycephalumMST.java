package edu.koritsas.slimemold.minimumtree;

import edu.koritsas.slimemold.AbstractPhysarumPolycephalum;
import edu.koritsas.slimemold.PhysarumPolycephalum;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

import java.util.List;

/**
 * Created by ilias on 5/10/2016.
 */
public abstract class PhysarumPolycephalumMST extends AbstractPhysarumPolycephalum {
    public PhysarumPolycephalumMST(Graph graph, Node sourceNode, List<Node> sinkNodes,double Io,double γ,int numberOfIterations){
        this.graph =graph;
        this.sourceNode=sourceNode;
        this.sinkNodesList=sinkNodes;
        this.Io=Io;
        this.γ=γ;
        this.numberOfIterations=numberOfIterations;

    }
    @Override
    public double calculateTubeConductivity(Edge e) {
        return 0;
    }
}
