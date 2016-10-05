package edu.koritsas.slimemold.minimumtree;

import edu.koritsas.slimemold.AbstractPhysarumPolycephalum;
import edu.koritsas.slimemold.PhysarumPolycephalum;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
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
        double D = conductivityMap.get(e);

        double Q = fluxMap.get(e);


        double fQ = Math.pow(Math.abs(Q), γ) / (1 + Math.pow(Math.abs(Q), γ));


        double newD = fQ - 0.1 * D;

       return newD;
    }

    @Override
    protected RealVector createConstantsVector() {
        List<Node> nodes = getAllNodes(graph);

        double[] constants =new double[nodes.size()];
        for (int i=0;i<nodes.size();i++){
            Node n=nodes.get(i);
            if (n.equals(sourceNode)) {
                constants[i] = Io;
            }else if(sinkNodesList.contains(n)){
                constants[i]=Io/(sinkNodesList.size());
            }else{
                constants[i]=0;
            }

        }


        return  new ArrayRealVector(constants);
    }
}
