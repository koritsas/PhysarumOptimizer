package edu.koritsas.slimemold;

import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.minimumtree.PhysarumPolycephalumMST;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

import java.util.List;

/**
 * Created by ilias on 5/10/2016.
 */
public class SlimeMold extends PhysarumPolycephalumMST {
    public SlimeMold(Graph graph, Node sourceNode, List<Node> sinkNodes, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNode, sinkNodes, Io, γ, numberOfIterations);
    }

    @Override
    public double getEdgeCost(Edge e) {
        SimpleFeature f = (SimpleFeature) e.getObject();
        Geometry geom = (Geometry) f.getDefaultGeometry();

        return geom.getLength();
    }
}
