import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.PhysarumPolycephalum;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

import java.util.List;

/**
 * Created by ilias on 29/9/2016.
 */
public class SlimeMold extends PhysarumPolycephalum {

    public SlimeMold(Graph graph, List<Node> sourceNodes, List<Node> sinkNodes, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNodes, sinkNodes, Io, γ, numberOfIterations);
    }

    @Override
    public double getEdgeCost(Edge e) {
        SimpleFeature f = (SimpleFeature) e.getObject();
        Geometry geometry = (Geometry) f.getDefaultGeometry();
        return geometry.getLength();
    }

    @Override
    public double getEdgeConstraint(Edge e) {
        return 0;
    }

    @Override
    public double getEuclideanDistance(Node n1, Node n2) {
        SimpleFeature f1 = (SimpleFeature) n1.getObject();
        SimpleFeature f2 = (SimpleFeature) n2.getObject();
        Geometry g1 = (Geometry) f1.getDefaultGeometry();
        Geometry g2 = (Geometry) f2.getDefaultGeometry();

        return g1.distance(g2);
    }
}
