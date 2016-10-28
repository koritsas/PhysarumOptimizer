package edu.koritsas.slimemold.shapefile;

import com.vividsolutions.jts.geom.Geometry;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.geotools.map.FeatureLayer;
import org.geotools.map.MapContent;
import org.geotools.styling.*;
import org.geotools.styling.Font;
import org.geotools.swing.JMapFrame;
import org.opengis.feature.simple.SimpleFeature;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ilias on 29/9/2016.
 */
public class GraphUtils {

    public static void visualizeGraph(Graph graph){
        List<Edge> edges =new ArrayList<Edge>(graph.getEdges());
        List<Node> nodes = new ArrayList<Node>(graph.getNodes());

     List<SimpleFeature> edgeFeatureList =new ArrayList<>();
        for (Edge e:edges){
            SimpleFeature f = (SimpleFeature) e.getObject();
            Geometry geometry = (Geometry) f.getDefaultGeometry();
            System.out.println(geometry);
            edgeFeatureList.add(f);
        }


    SimpleFeatureCollection edgeFeatureCollection = new ListFeatureCollection(edgeFeatureList.get(0).getFeatureType(),edgeFeatureList);

        Style edgeStyle =SLD.createLineStyle(Color.green,0.5f,"SHAPE_Leng",null);
        FeatureLayer edgeLayer = new FeatureLayer(edgeFeatureCollection,edgeStyle);

   List<SimpleFeature> nodeFeature =new ArrayList<>();
        for (Node n:nodes){
            try {
                SimpleFeature feature = (SimpleFeature) n.getObject();
                nodeFeature.add(feature);
            }catch (ClassCastException e){
                System.out.println("For the last time! It's not  feature, it's a point!!");
            }
        }
      SimpleFeatureCollection nodeFeatureCollection =new ListFeatureCollection(nodeFeature.get(0).getFeatureType(),nodeFeature);
        Style nodeStyle =SLD.createPointStyle("Circle",Color.black,Color.gray,0f,5f,"OBJECTID", null);

        FeatureLayer nodeLayer = new FeatureLayer(nodeFeatureCollection,nodeStyle);

        MapContent mapContent = new MapContent();
        mapContent.addLayer(edgeLayer);
       mapContent.addLayer(nodeLayer);

        JMapFrame mapFrame=new JMapFrame();
        mapFrame.showMap(mapContent);



    }




}
