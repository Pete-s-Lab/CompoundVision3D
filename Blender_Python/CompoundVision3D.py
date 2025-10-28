bl_info = {
    "name": "CompoundVision3D",
    "blender": (3, 0, 0),
    "category": "View3D",
    "author": "Peter T. Rühr",
    "version": (0, 0, 9023),
    "description": "CompoundVision3D Blender Plugin. Part of the compound eye property extraction workflow published by Rühr, Pande & Blanke: xxx.",
}

import bpy
import csv
import os.path as p
import re

        
# ------------------ SCRIPT 1 ------------------
class SCRIPT_1_OT_Operator(bpy.types.Operator):
    bl_idname = "script.execute_1"
    bl_label = "Step 1: Rescale mesh"
    
    def execute(self, context):
        bpy.ops.transform.resize(value=(0.001, 0.001, 0.001), orient_type='GLOBAL')
        bpy.ops.object.editmode_toggle()
        bpy.ops.mesh.normals_make_consistent(inside=False)
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.flip_normals()
        bpy.ops.mesh.select_all(action='DESELECT')
        for a in bpy.context.screen.areas:
            if a.type == 'VIEW_3D':
                for s in a.spaces:
                    if s.type == 'VIEW_3D':
                        s.clip_end = 10000
                        s.overlay.show_stats = True
        self.report({'INFO'}, "Script 1 executed!")
        return {'FINISHED'}

# ------------------ SCRIPT 2 ------------------
class SCRIPT_2_OT_Operator(bpy.types.Operator):
    bl_idname = "script.execute_2"
    bl_label = "Step 2: Add modifiers before export"
    
    def execute(self, context):
        bpy.ops.object.modifier_add(type='DECIMATE')
        bpy.context.object.modifiers["Decimate"].ratio = 0.5  # Example ratio
        bpy.ops.object.modifier_add(type='SMOOTH')
        bpy.context.object.modifiers["Smooth"].iterations = 10
        bpy.context.object.modifiers["Smooth"].factor = 0.5
        bpy.ops.object.editmode_toggle()
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.flip_normals()
        bpy.ops.object.editmode_toggle()
        self.report({'INFO'}, "Script 2 executed!")
        return {'FINISHED'}

# ------------------ SCRIPT 3 ------------------
class SCRIPT_3_OT_Operator(bpy.types.Operator):
    bl_idname = "script.execute_3"
    bl_label = "Run Script 3"
    
    def execute(self, context):
        # Step 3 Settings:
        diameter = 2.5/4 # 25 2.5
        segments = 16 # 4 16
        ring_count = 8 # 3 8 

        camera_type = 'PERSP' # 'ORTHO' 'PERSP'
        
        curr_filename_raw = bpy.path.basename(bpy.context.blend_data.filepath)
        curr_filename = re.sub(".blend$", "", curr_filename_raw)
        curr_filepath = bpy.path.abspath("//")
        candidates_filepath = re.sub("1_pre_STLs.+$", "6_facet_candidates", curr_filepath)
        
        csv_file_path = p.join(candidates_filepath, curr_filename+'_facet_candidates.csv')
        file = csv.reader(open(csv_file_path, newline=''), delimiter=',')
        radius = 0.01
        
        
        # Function to create a sphere template
        def create_sphere_template(radius=0.1):
            mesh = bpy.data.meshes.new(name="TemplateSphereMesh")
            sphere = bpy.data.objects.new(name="TemplateSphere", object_data=mesh)
            bpy.context.collection.objects.link(sphere)

            # Create a UV sphere mesh
            bpy.ops.object.select_all(action='DESELECT')
            bpy.context.view_layer.objects.active = sphere
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.primitive_uv_sphere_add(segments=16, ring_count=8, radius=radius)
            bpy.ops.object.mode_set(mode='OBJECT')
            
            return sphere

        # Function to add spheres for points from a CSV
        def add_spheres_from_csv(filepath, radius=10):
            # Create the template sphere
            template_sphere = create_sphere_template(radius)

            # Read points from the CSV file
            with open(filepath, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    # Extract coordinates and ID
                    x, y, z = float(row['x'])/1000, float(row['y'])/1000, float(row['z'])/1000
                    point_id = row['ID']
                    
                    # Create an instance of the template sphere
                    sphere_instance = template_sphere.copy()
                    sphere_instance.location = (x, y, z)
                    sphere_instance.name = f"{point_id}"  # Name the sphere based on the ID
                    bpy.context.collection.objects.link(sphere_instance)
            
            # Delete the original template sphere to save memory
            bpy.data.objects.remove(template_sphere)

        # Add spheres for the point cloud in the CSV
        add_spheres_from_csv(csv_file_path, radius=diameter/2/100)
        
        
        # for idx, row in enumerate(file):
            # if idx > 0:
                # x = float(row[1])
                # y = float(row[2])
                # z = float(row[3])
                # bpy.ops.object.empty_add(type='SPHERE', radius=radius, align='WORLD', location=(x/1000, y/1000, z/1000), scale=(1, 1, 1))
        self.report({'INFO'}, "Script 3 executed!")
        return {'FINISHED'}

# ------------------ SCRIPT 4 ------------------
class SCRIPT_4_OT_Operator(bpy.types.Operator):
    bl_idname = "script.execute_4"
    bl_label = "Run Script 4"
    
    def execute(self, context):
        curr_filename_raw = bpy.path.basename(bpy.context.blend_data.filepath)
        curr_filename = re.sub(".blend$", "", curr_filename_raw)
        curr_filepath = bpy.path.abspath("//")
        positions_filepath = re.sub("1_pre_STLs.+$", "7_facet_positions", curr_filepath)
        
        coordinate_multiplicator = 1000
        
        output_path = p.join(positions_filepath, curr_filename+'_facet_positions.csv')
        
        selected_objects = bpy.context.selected_objects
        with open(output_path, "w") as output:
            output.write("ID,x,y,z\n")
            for obj in selected_objects:
                string_to_write = (f"{obj.name}", str(obj.location[0]*coordinate_multiplicator), str(obj.location[1]*coordinate_multiplicator), str(obj.location[2]*coordinate_multiplicator), "\n")
                output.write(','.join(string_to_write))
        self.report({'INFO'}, f"Data exported to: {output_path}")
        return {'FINISHED'}

# ------------------ GUI PANEL ------------------
class MYPLUGIN_PT_Panel(bpy.types.Panel):
    bl_label = "CompoundVision3D Panel"
    bl_idname = "MYPLUGIN_PT_Panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "CompoundVision3D"
    
    def draw(self, context):
        layout = self.layout
        layout.label(text="Execute Steps:")
        layout.operator("script.execute_1", text="1: Rescale mesh")
        layout.operator("script.execute_2", text="2: Add modifiers")
        layout.operator("script.execute_3", text="3: Import candidates")
        layout.operator("script.execute_4", text="4: Export positions")

# ------------------ REGISTER ADD-ON ------------------
classes = [
    SCRIPT_1_OT_Operator,
    SCRIPT_2_OT_Operator,
    SCRIPT_3_OT_Operator,
    SCRIPT_4_OT_Operator,
    MYPLUGIN_PT_Panel,
]

def register():
    for cls in classes:
        bpy.utils.register_class(cls)

def unregister():
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

if __name__ == "__main__":
    register()
