from flask import Flask
from flask_cors import CORS


def create_app():
    app = Flask(__name__)
    app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024  # 50 MB upload limit
    CORS(app)

    from routes.process import bp as process_bp
    from routes.meta import bp as meta_bp

    app.register_blueprint(process_bp)
    app.register_blueprint(meta_bp)

    return app


if __name__ == "__main__":
    create_app().run(debug=True, port=5000)
